use super::*;
use crate::leb::LebTier;
use crate::positions::{EntrySnapshot, StateProvider};
use crate::types::TwoPartEpoch;

/// Dummy provider that never serves anything: `has_body` is false so
/// the sealed gate refuses before any state read; every read is a
/// typed miss. Enough for the binding-level refusal tests below,
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

#[test]
fn jd_tt_adds_deltat_and_refuses_nonfinite() {
    let provider = Unserved;
    let nut = |_: f64| Ok((0.0, 0.0));
    let serving = SealedServing::new(&provider, &nut);
    let dt = |jd: f64| Ok(jd * 0.0 + 0.001);
    let bound = SealedLongitudeServing::new(
        &serving,
        &dt,
        crate::constants::SUN,
        crate::constants::FLG_SWIEPH,
    );
    assert!((bound.jd_tt_for_ut(2460000.0).unwrap() - 2460000.001).abs() < 1e-12);
    assert!(matches!(
        bound.jd_tt_for_ut(f64::NAN),
        Err(Error::NonFiniteArgument { .. })
    ));
    let bad_dt = |_: f64| Ok(f64::INFINITY);
    let bad = SealedLongitudeServing::new(&serving, &bad_dt, crate::constants::SUN, 2);
    assert!(matches!(
        bad.jd_tt_for_ut(2460000.0),
        Err(Error::NonFiniteArgument { .. })
    ));
    let _ = zero_dt(0.0);
}

#[test]
fn lon_speed_propagates_sealed_miss_fail_closed() {
    let provider = Unserved;
    let nut = |_: f64| Ok((0.0, 0.0));
    let serving = SealedServing::new(&provider, &nut);
    let dt = |_: f64| Ok(0.0008);
    let bound = SealedLongitudeServing::new(
        &serving,
        &dt,
        crate::constants::SUN,
        crate::constants::FLG_SWIEPH,
    );
    // Unserved channel: the sealed gate refuses with a typed miss;
    // the walks propagate it fail-closed, never as a fabricated
    // crossing (range misses become exhaustion inside the walks).
    let err = bound.lon_speed_at_ut(2460000.0).unwrap_err();
    assert!(matches!(err, Error::LebBodyNotFound { .. }));
}

#[test]
fn station_search_refuses_sun_without_provider_contact() {
    let provider = Unserved;
    let nut = |_: f64| Ok((0.0, 0.0));
    let serving = SealedServing::new(&provider, &nut);
    let dt = |_: f64| Ok(0.0008);
    let bound = SealedLongitudeServing::new(
        &serving,
        &dt,
        crate::constants::MARS,
        crate::constants::FLG_SWIEPH,
    );
    // Sun/Moon never station: refused before any
    // state read, even though the bound body is a planet.
    assert!(matches!(
        find_station_ut_serving(
            &bound,
            crate::constants::SUN,
            2460000.0,
            StationRequest::Any
        ),
        Err(Error::BodyHasNoStations { body: 0 })
    ));
    assert!(matches!(
        next_retrograde_ut_serving(&bound, crate::constants::MOON, 2460000.0),
        Err(Error::BodyHasNoStations { body: 1 })
    ));
}

#[test]
fn station_info_servings_fold_sun_moon_and_misses() {
    // Release `is_retrograde` / `get_station_type` read `false` /
    // `"direct"` for Sun/Moon without consulting the ephemeris, and
    // fold every provider failure the same way (never a raised
    // error). `station_info` / `velocity_at_station` stay fail-closed.
    let provider = Unserved;
    let nut = |_: f64| Ok((0.0, 0.0));
    let serving = SealedServing::new(&provider, &nut);
    let dt = |_: f64| Ok(0.0008);
    let bound = SealedLongitudeServing::new(
        &serving,
        &dt,
        crate::constants::MARS,
        crate::constants::FLG_SWIEPH,
    );
    assert!(!is_retrograde_ut_serving(&bound, crate::constants::SUN, 2460000.0).unwrap());
    assert!(!is_retrograde_ut_serving(&bound, crate::constants::MOON, 2460000.0).unwrap());
    // Unserved channel: the miss folds to `false`, never an error.
    assert!(!is_retrograde_ut_serving(&bound, crate::constants::MARS, 2460000.0).unwrap());
    assert_eq!(
        station_type_ut_serving(&bound, crate::constants::SUN, 2460000.0).unwrap(),
        StationMotionKind::Direct
    );
    assert_eq!(
        station_type_ut_serving(&bound, crate::constants::MOON, 2460000.0).unwrap(),
        StationMotionKind::Direct
    );
    assert_eq!(
        station_type_ut_serving(&bound, crate::constants::MARS, 2460000.0).unwrap(),
        StationMotionKind::Direct
    );
    assert_eq!(StationMotionKind::Direct.code(), 0);
    assert_eq!(StationMotionKind::Retrograde.code(), 1);
    assert_eq!(StationMotionKind::StationaryRetrograde.code(), 2);
    assert_eq!(StationMotionKind::StationaryDirect.code(), 3);
    assert_eq!(
        station_type_ut_serving(&bound, crate::constants::MARS, f64::NAN).unwrap(),
        StationMotionKind::Direct
    );
    // `station_info` refuses Sun/Moon before any state read and
    // propagates the sealed miss fail-closed (never a fabricated
    // station); the velocity arm refuses a non-positive step
    // before any state read.
    assert!(matches!(
        station_info_ut_serving(&bound, crate::constants::SUN, 2460000.0),
        Err(Error::BodyHasNoStations { body: 0 })
    ));
    assert!(matches!(
        station_info_ut_serving(&bound, crate::constants::MOON, 2460000.0),
        Err(Error::BodyHasNoStations { body: 1 })
    ));
    assert!(matches!(
        station_info_ut_serving(&bound, crate::constants::MARS, 2460000.0),
        Err(Error::LebBodyNotFound { .. })
    ));
    assert!(matches!(
        station_info_ut_serving(&bound, crate::constants::MARS, f64::NAN),
        Err(Error::NonFiniteSearchValue { .. })
    ));
    assert!(matches!(
        velocity_at_station_ut_serving(&bound, 2460000.0, 1.0),
        Err(Error::LebBodyNotFound { .. })
    ));
    assert!(matches!(
        velocity_at_station_ut_serving(&bound, 2460000.0, 0.0),
        Err(Error::InvalidSearchStep { .. })
    ));
    assert_eq!(
        StationMotionKind::StationaryRetrograde.as_str(),
        "stationary_retrograde"
    );
}

#[test]
fn crossing_wrappers_refuse_nonfinite_start() {
    let provider = Unserved;
    let nut = |_: f64| Ok((0.0, 0.0));
    let serving = SealedServing::new(&provider, &nut);
    let dt = |_: f64| Ok(0.0008);
    let bound = SealedLongitudeServing::new(
        &serving,
        &dt,
        crate::constants::SUN,
        crate::constants::FLG_SWIEPH,
    );
    assert!(matches!(
        solcross_ut_serving(&bound, 0.0, f64::NAN, SearchDirection::Forward),
        Err(Error::NonFiniteSearchValue { .. })
    ));
    assert!(matches!(
        mooncross_ut_serving(&bound, 0.0, f64::INFINITY, SearchDirection::Forward),
        Err(Error::NonFiniteSearchValue { .. })
    ));
    assert!(matches!(
        cross_ut_serving(&bound, crate::constants::MARS, 0.0, f64::NAN),
        Err(Error::NonFiniteSearchValue { .. })
    ));
    assert!(matches!(
        mooncross_node_ut_serving(&bound, f64::NAN, SearchDirection::Forward),
        Err(Error::NonFiniteSearchValue { .. })
    ));
}

#[test]
fn angular_servings_refuse_format_flags_fix_c34() {
    // FIX-C34: FLG_RADIANS/FLG_XYZ silently corrupt the
    // angular-word reading (radian/Cartesian words vs degree
    // targets), so every longitude/transit entry refuses them
    // before any provider contact.
    let provider = Unserved;
    let nut = |_: f64| Ok((0.0, 0.0));
    let serving = SealedServing::new(&provider, &nut);
    let dt = |_: f64| Ok(0.0008);
    let rad = SealedLongitudeServing::new(
        &serving,
        &dt,
        crate::constants::SUN,
        crate::constants::FLG_RADIANS,
    );
    assert!(matches!(
        solcross_ut_serving(&rad, 0.0, 2460000.0, SearchDirection::Forward),
        Err(Error::UnsupportedFlagCombination { .. })
    ));
    assert!(matches!(
        mooncross_node_ut_serving(&rad, 2460000.0, SearchDirection::Forward),
        Err(Error::UnsupportedFlagCombination { .. })
    ));
    let xyz = SealedLongitudeServing::new(
        &serving,
        &dt,
        crate::constants::MARS,
        crate::constants::FLG_XYZ,
    );
    assert!(matches!(
        cross_ut_serving(&xyz, crate::constants::MARS, 0.0, 2460000.0),
        Err(Error::UnsupportedFlagCombination { .. })
    ));
    assert!(matches!(
        find_station_ut_serving(&xyz, crate::constants::MARS, 2460000.0, StationRequest::Any),
        Err(Error::UnsupportedFlagCombination { .. })
    ));
    assert!(matches!(
        next_retrograde_ut_serving(&xyz, crate::constants::MARS, 2460000.0),
        Err(Error::UnsupportedFlagCombination { .. })
    ));
}

#[test]
fn rise_set_serving_refuses_format_flags_fix_c37() {
    // FIX-C37: the rise/set altitude engine reads angular RA/Dec
    // words from the transit serving, so radians/XYZ servings are
    // refused on this path too (the FIX-C34 guard covers only the
    // culmination entry).
    let provider = Unserved;
    let nut = |_: f64| Ok((0.0, 0.0));
    let serving = SealedServing::new(&provider, &nut);
    let dt = |_: f64| Ok(0.0008);
    let rad = SealedTransitServing::try_new(
        &serving,
        &dt,
        crate::constants::SUN,
        crate::constants::FLG_RADIANS,
    )
    .expect("offsets resolve");
    let req = RiseSetRequest {
        jd_start_ut: 2460000.0,
        lon_deg_east: 12.5,
        lat_deg: 41.9,
        height_m: 0.0,
        atpress_mbar: 1013.25,
        attemp_c: 15.0,
        horhgt_deg: 0.0,
        rsmi: crate::constants::CALC_RISE,
    };
    assert!(matches!(
        rise_set_ut_serving(&rad, &req),
        Err(Error::UnsupportedFlagCombination { .. })
    ));
}

// -- slice 1: sealed transit serving tests --------------------

fn transit_bound<'a>(
    serving: &'a SealedServing<'a, Unserved>,
    dt: &'a dyn Fn(f64) -> Result<f64>,
) -> SealedTransitServing<'a, Unserved> {
    SealedTransitServing::try_new(
        serving,
        dt,
        crate::constants::SUN,
        crate::constants::FLG_SWIEPH,
    )
    .expect("offsets resolve")
}

#[test]
fn transit_serving_resolves_offsets_and_adds_deltat() {
    let provider = Unserved;
    let nut = |_: f64| Ok((0.0, 0.0));
    let serving = SealedServing::new(&provider, &nut);
    let dt = |_: f64| Ok(0.0008);
    let bound = transit_bound(&serving, &dt);
    // Boundary offsets match the direct composition over the
    // same provider words (wiring pin, not a model claim).
    use crate::time::sidereal_continuity_offsets;
    let (lo, hi) = sidereal_continuity_offsets(0.0008, 0.0008).unwrap();
    assert_eq!(bound.offset_lo_deg, lo);
    assert_eq!(bound.offset_hi_deg, hi);
    assert!((bound.jd_tt_for_ut(2460000.0).unwrap() - 2460000.0008).abs() < 1e-12);
    assert!(matches!(
        bound.jd_tt_for_ut(f64::NAN),
        Err(Error::NonFiniteArgument { .. })
    ));
}

#[test]
fn transit_serving_refuses_bad_deltat_words() {
    let provider = Unserved;
    let nut = |_: f64| Ok((0.0, 0.0));
    let serving = SealedServing::new(&provider, &nut);
    let bad_dt = |_: f64| Ok(f64::INFINITY);
    assert!(matches!(
        SealedTransitServing::try_new(
            &serving,
            &bad_dt,
            crate::constants::SUN,
            crate::constants::FLG_SWIEPH,
        ),
        Err(Error::NonFiniteArgument { .. })
    ));
    let dt = |_: f64| Ok(0.0008);
    let bound = transit_bound(&serving, &dt);
    let bad = |_: f64| Ok(f64::NAN);
    let bad_bound = SealedTransitServing::try_new(
        &serving,
        &bad,
        crate::constants::SUN,
        crate::constants::FLG_SWIEPH,
    );
    assert!(bad_bound.is_err());
    assert!(matches!(
        bound.gast_hours_at_ut(f64::NAN),
        Err(Error::NonFiniteArgument { .. })
    ));
}

#[test]
fn transit_serving_threads_nutation_and_obliquity_into_gast() {
    use crate::time::{SiderealDeltaT, greenwich_apparent_sidereal_time_hours};
    let provider = Unserved;
    // 10 mas dpsi: a non-zero word so a rad/deg mixup or a dropped
    // term moves the answer far outside exact equality.
    let nut = |_: f64| Ok(((10.0_f64 / 3_600_000.0).to_radians(), 0.0));
    let serving = SealedServing::new(&provider, &nut);
    let dt = |_: f64| Ok(0.0008);
    let bound = transit_bound(&serving, &dt);
    let jd_ut = 2460000.0;
    let jd_tt = jd_ut + 0.0008;
    let (dpsi_rad, deps_rad) = nut(jd_tt).unwrap();
    let eps_true =
        crate::positions::vondrak_mean_obliquity_deg(jd_tt).unwrap() + deps_rad.to_degrees();
    let sdt = SiderealDeltaT::try_new(0.0008, bound.offset_lo_deg, bound.offset_hi_deg).unwrap();
    let direct =
        greenwich_apparent_sidereal_time_hours(jd_ut, eps_true, dpsi_rad.to_degrees(), &sdt)
            .unwrap();
    assert_eq!(bound.gast_hours_at_ut(jd_ut).unwrap(), direct);
}

#[test]
fn transit_serving_propagates_sealed_miss_fail_closed() {
    let provider = Unserved;
    let nut = |_: f64| Ok((0.0, 0.0));
    let serving = SealedServing::new(&provider, &nut);
    let dt = |_: f64| Ok(0.0008);
    let bound = transit_bound(&serving, &dt);
    // Unserved channel: the sealed gate refuses with a typed miss;
    // the transit search propagates it fail-closed, never as a
    // fabricated transit.
    let err = bound.ra_dec_at_ut(2460000.0).unwrap_err();
    assert!(matches!(err, Error::LebBodyNotFound { .. }));
    let err = bound.transit_at_ut(2460000.0, 0.0, false).unwrap_err();
    assert!(matches!(err, Error::LebBodyNotFound { .. }));
}

#[test]
fn transit_wrappers_refuse_nonfinite_start_and_lon() {
    let provider = Unserved;
    let nut = |_: f64| Ok((0.0, 0.0));
    let serving = SealedServing::new(&provider, &nut);
    let dt = |_: f64| Ok(0.0008);
    let bound = transit_bound(&serving, &dt);
    assert!(matches!(
        transit_ut_serving(&bound, f64::NAN, 0.0, false),
        Err(Error::NonFiniteSearchValue { .. })
    ));
    assert!(matches!(
        transit_ut_serving(&bound, 2460000.0, f64::INFINITY, true),
        Err(Error::InvalidObserver { .. })
    ));
}

// -- slice 2: sealed rise/set serving tests -------------------

#[test]
fn cotrans_ecliptic_equatorial_round_trip() {
    // Release `utils.cotrans` docstring pins (negative eps maps
    // ecliptic to equatorial and back).
    let (ra, dec, dist) = cotrans(90.0, 0.0, 1.0, -23.4);
    assert!((ra - 90.0).abs() < 1e-12);
    assert!((dec - 23.4).abs() < 1e-12);
    assert_eq!(dist, 1.0);
    let (lon, lat, dist) = cotrans(90.0, 23.4, 1.0, 23.4);
    assert!((lon - 90.0).abs() < 1e-12);
    assert!(lat.abs() < 1e-12);
    assert_eq!(dist, 1.0);
    // RA normalization: a negative branch lands in [0, 360).
    let (ra_neg, _, _) = cotrans(270.0, 0.0, 1.0, -23.4);
    assert!((0.0..360.0).contains(&ra_neg));
}

#[test]
fn cotrans_degree_and_radian_entries_preserve_nonlongitude_words() {
    for lon in [0.0_f64, -0.0] {
        for lat in [0.0_f64, -0.0] {
            for eps in [0.0_f64, -0.0, 360.0, -360.0] {
                for distance in [0.0_f64, -0.0, 1.0, -1.25] {
                    let degrees = cotrans(lon, lat, distance, eps);
                    let radians = crate::events::cotrans_radians(
                        lon.to_radians(),
                        lat.to_radians(),
                        (-eps).to_radians(),
                        distance,
                    );
                    for result in [degrees, radians] {
                        assert_eq!(result.0.to_bits(), 0.0_f64.to_bits());
                        assert_eq!(result.1, 0.0);
                        assert_eq!(result.2.to_bits(), distance.to_bits());
                    }
                    assert_eq!(degrees.1.to_bits(), radians.1.to_bits());
                }
            }
        }
    }
    assert_eq!(
        cotrans(0.0, -0.0, 1.0, 0.0).1.to_bits(),
        (-0.0_f64).to_bits()
    );
    assert_eq!(cotrans(-0.0, -0.0, 1.0, 0.0).1.to_bits(), 0.0_f64.to_bits());
    for lon in [1e-14_f64, -1e-14, 360.0, -360.0, 720.0, -720.0] {
        assert_ne!(cotrans(lon, 0.0, 1.0, 0.0).0, 0.0);
        assert_ne!(
            crate::events::cotrans_radians(lon.to_radians(), 0.0, -0.0, 1.0).0,
            0.0
        );
    }
}

#[test]
fn body_disc_radii_match_release_expressions() {
    assert_eq!(
        SUN_RADIUS_AU.to_bits(),
        (696000.0f64 / 149597870.700f64).to_bits()
    );
    assert_eq!(
        MOON_RADIUS_AU.to_bits(),
        (1738.15f64 / 149597870.700f64).to_bits()
    );
    assert_eq!(
        body_disc_radius_au(crate::constants::SUN).to_bits(),
        SUN_RADIUS_AU.to_bits()
    );
    assert_eq!(
        body_disc_radius_au(crate::constants::MOON).to_bits(),
        MOON_RADIUS_AU.to_bits()
    );
    assert_eq!(
        body_disc_radius_au(crate::constants::MARS).to_bits(),
        (3389.5f64 / 149597870.700f64).to_bits()
    );
    // No finite disc on this path: nodes, fictitious points, unknown ids.
    assert_eq!(body_disc_radius_au(crate::constants::ECL_NUT), 0.0);
    assert_eq!(body_disc_radius_au(99999), 0.0);
}

fn rise_request(rsmi: i32) -> RiseSetRequest {
    RiseSetRequest {
        jd_start_ut: 2460000.0,
        lon_deg_east: 12.5,
        lat_deg: 41.9,
        height_m: 0.0,
        atpress_mbar: 1013.25,
        attemp_c: 15.0,
        horhgt_deg: 0.0,
        rsmi,
    }
}

#[test]
fn rise_set_refuses_degenerate_bodies() {
    use crate::constants::{EARTH, ECL_NUT, FLG_SWIEPH};
    let provider = Unserved;
    let nut = |_: f64| Ok((0.0, 0.0));
    let serving = SealedServing::new(&provider, &nut);
    let dt = |_: f64| Ok(0.0008);
    for body in [EARTH, ECL_NUT] {
        let bound = SealedTransitServing::try_new(&serving, &dt, body, FLG_SWIEPH)
            .expect("offsets resolve");
        assert!(matches!(
            rise_set_ut_serving(&bound, &rise_request(crate::constants::CALC_RISE)),
            Err(Error::InvalidBody { body_id, .. }) if body_id == body
        ));
    }
}

#[test]
fn rise_set_refuses_transit_kinds_and_nonfinite_inputs() {
    use crate::constants::{CALC_ITRANSIT, CALC_MTRANSIT, CALC_RISE};
    let provider = Unserved;
    let nut = |_: f64| Ok((0.0, 0.0));
    let serving = SealedServing::new(&provider, &nut);
    let dt = |_: f64| Ok(0.0008);
    let bound = transit_bound(&serving, &dt);
    for kind in [CALC_MTRANSIT, CALC_ITRANSIT] {
        assert!(matches!(
            rise_set_ut_serving(&bound, &rise_request(kind)),
            Err(Error::InvalidSearchStep { .. })
        ));
    }
    let mut bad_start = rise_request(CALC_RISE);
    bad_start.jd_start_ut = f64::NAN;
    assert!(matches!(
        rise_set_ut_serving(&bound, &bad_start),
        Err(Error::NonFiniteSearchValue { .. })
    ));
    let mut bad_lon = rise_request(CALC_RISE);
    bad_lon.lon_deg_east = f64::INFINITY;
    assert!(matches!(
        rise_set_ut_serving(&bound, &bad_lon),
        Err(Error::InvalidObserver { .. })
    ));
    assert!(matches!(
        bound.last_hours_at_ut(f64::NAN, 0.0),
        Err(Error::NonFiniteArgument { .. })
    ));
    assert!(matches!(
        bound.last_hours_at_ut(2460000.0, f64::NAN),
        Err(Error::InvalidObserver { .. })
    ));
}

#[test]
fn rise_set_propagates_sealed_miss_fail_closed() {
    use crate::constants::CALC_SET;
    let provider = Unserved;
    let nut = |_: f64| Ok((0.0, 0.0));
    let serving = SealedServing::new(&provider, &nut);
    let dt = |_: f64| Ok(0.0008);
    let bound = transit_bound(&serving, &dt);
    // Unserved channel: the sealed gate refuses with a typed miss;
    // the rise/set search propagates it fail-closed, never as a
    // fabricated event or a silent circumpolar verdict.
    let err = rise_set_ut_serving(&bound, &rise_request(CALC_SET)).unwrap_err();
    assert!(matches!(err, Error::LebBodyNotFound { .. }));
}

// -- slice 3: sealed star + Hindu tests ----------------------

fn star_request(name_rsmi: i32) -> StarRiseSetRequest {
    StarRiseSetRequest {
        jd_start_ut: 2460000.0,
        lon_deg_east: 12.5,
        lat_deg: 41.9,
        height_m: 0.0,
        atpress_mbar: 1013.25,
        attemp_c: 15.0,
        horhgt_deg: 0.0,
        rsmi: name_rsmi,
    }
}

fn star_bound<'a>(
    star_serving: &'a crate::stars::SealedStarServing<'a, Unserved>,
    dt: &'a dyn Fn(f64) -> Result<f64>,
    nut: &'a dyn Fn(f64) -> Result<(f64, f64)>,
) -> SealedStarEventServing<'a, Unserved> {
    SealedStarEventServing::try_new(star_serving, "Sirius", dt, nut).expect("Sirius resolves")
}

#[test]
fn star_serving_resolves_v1_names_and_refuses_unknown() {
    use crate::sidereal::SiderealSelection;
    let provider = Unserved;
    let nut = |_: f64| Ok((0.0, 0.0));
    let dt = |_: f64| Ok(0.0008);
    let star_serving =
        crate::stars::SealedStarServing::new(&provider, &nut, &dt, SiderealSelection::unset());
    // v1 grammar: a canonical bright-star name resolves with the
    // release `"Name,nomenclature"` canonical form.
    let bound = SealedStarEventServing::try_new(&star_serving, "Sirius", &dt, &nut).unwrap();
    assert_eq!(bound.star().canonical, "Sirius,alCMa");
    // Unknown names refuse StarNotFound (never a silent substitution);
    // the release `int` digit edge maps to the typed refusal.
    assert!(matches!(
        SealedStarEventServing::try_new(&star_serving, "No Such Star xyz", &dt, &nut),
        Err(Error::StarNotFound { .. })
    ));
    // Non-finite Delta-T words refuse at construction.
    let bad_dt = |_: f64| Ok(f64::NAN);
    assert!(matches!(
        SealedStarEventServing::try_new(&star_serving, "Sirius", &bad_dt, &nut),
        Err(Error::NonFiniteArgument { .. })
    ));
}

#[test]
fn star_disc_is_always_center() {
    use crate::constants::{
        BIT_DISC_BOTTOM, BIT_DISC_CENTER, BIT_FIXED_DISC_SIZE, BIT_NO_REFRACTION, CALC_RISE,
    };
    // Disc modifier bits are accepted and ignored for fixed stars:
    // the release always uses the disc center (`disc_sign = 0`).
    for bits in [
        0,
        BIT_DISC_CENTER,
        BIT_DISC_BOTTOM,
        BIT_FIXED_DISC_SIZE,
        BIT_NO_REFRACTION,
        BIT_DISC_BOTTOM | BIT_FIXED_DISC_SIZE,
    ] {
        assert_eq!(
            select_disc_point(CALC_RISE | bits, true),
            DiscPoint::Center,
            "bits {bits:#x}"
        );
    }
}

#[test]
fn star_rise_set_validates_before_any_ephemeris_read() {
    use crate::constants::{CALC_ITRANSIT, CALC_MTRANSIT, CALC_RISE};
    use crate::sidereal::SiderealSelection;
    let provider = Unserved;
    let nut = |_: f64| Ok((0.0, 0.0));
    let dt = |_: f64| Ok(0.0008);
    let star_serving =
        crate::stars::SealedStarServing::new(&provider, &nut, &dt, SiderealSelection::unset());
    let bound = star_bound(&star_serving, &dt, &nut);
    // Transit kinds belong to the star transit entry.
    for kind in [CALC_MTRANSIT, CALC_ITRANSIT] {
        assert!(matches!(
            star_rise_set_ut_serving(&bound, &star_request(kind)),
            Err(Error::InvalidSearchStep { .. })
        ));
    }
    let mut bad_start = star_request(CALC_RISE);
    bad_start.jd_start_ut = f64::NAN;
    assert!(matches!(
        star_rise_set_ut_serving(&bound, &bad_start),
        Err(Error::NonFiniteSearchValue { .. })
    ));
    let mut bad_lon = star_request(CALC_RISE);
    bad_lon.lon_deg_east = f64::INFINITY;
    assert!(matches!(
        star_rise_set_ut_serving(&bound, &bad_lon),
        Err(Error::InvalidObserver { .. })
    ));
    // The unserved provider's typed miss propagates fail-closed past
    // validation (never a fabricated event or silent circumpolar).
    let err = star_rise_set_ut_serving(&bound, &star_request(CALC_RISE)).unwrap_err();
    assert!(matches!(err, Error::LebBodyNotFound { .. }));
}

#[test]
fn star_events_ignore_sidereal_selection() {
    use crate::constants::{CALC_RISE, SIDM_TRUE_CITRA};
    use crate::sidereal::SiderealSelection;
    // STAR-04 event-integration disposition: the release star-event
    // acquisition never forwards sidereal flags (tropical `FLG_SPEED` /
    // `FLG_EQUATORIAL` only; the pin `eclipse.py` references no sidereal
    // mode, and occult servings refuse `FLG_SIDEREAL` at construction),
    // so the bound sidereal selection must not move star-event outcomes
    // and dynamic live anchors need no event-caller wiring. A dynamic
    // selection that serves live anchor places must behave exactly like
    // the unset default on this path; any future anchor consultation in
    // the tropical acquisition would surface here as a divergence from
    // the unserved-provider miss.
    let provider = Unserved;
    let nut = |_: f64| Ok((0.0, 0.0));
    let dt = |_: f64| Ok(0.0008);
    for (label, selection) in [
        ("unset", SiderealSelection::unset()),
        ("dynamic", {
            let mut sel = SiderealSelection::unset();
            let _warnings = sel.set_sidereal_mode(SIDM_TRUE_CITRA, 0.0, 0.0);
            sel
        }),
    ] {
        let star_serving = crate::stars::SealedStarServing::new(&provider, &nut, &dt, selection);
        let bound = star_bound(&star_serving, &dt, &nut);
        let err = star_rise_set_ut_serving(&bound, &star_request(CALC_RISE)).unwrap_err();
        assert!(
            matches!(err, Error::LebBodyNotFound { .. }),
            "{label}: expected the unserved-provider miss, got {err:?}"
        );
    }
}

#[test]
fn hindu_branch_propagates_sealed_miss_fail_closed() {
    use crate::constants::{BIT_GEOCTR_NO_ECL_LAT, CALC_RISE, FLG_SWIEPH};
    let provider = Unserved;
    let nut = |_: f64| Ok((0.0, 0.0));
    let serving = SealedServing::new(&provider, &nut);
    let dt = |_: f64| Ok(0.0008);
    let bound = SealedTransitServing::try_new(&serving, &dt, crate::constants::MARS, FLG_SWIEPH)
        .expect("offsets resolve");
    // The Hindu geocentric fetch hits the same sealed gate: an
    // unserved channel refuses with the typed miss on both the
    // ordinary and the Hindu path.
    let err =
        rise_set_ut_serving(&bound, &rise_request(CALC_RISE | BIT_GEOCTR_NO_ECL_LAT)).unwrap_err();
    assert!(matches!(err, Error::LebBodyNotFound { .. }));
}

#[test]
fn star_gast_threads_nutation_and_obliquity() {
    use crate::sidereal::SiderealSelection;
    use crate::time::{SiderealDeltaT, greenwich_apparent_sidereal_time_hours};
    let provider = Unserved;
    let nut = |_: f64| Ok(((10.0_f64 / 3_600_000.0).to_radians(), 0.0));
    let dt = |_: f64| Ok(0.0008);
    let star_serving =
        crate::stars::SealedStarServing::new(&provider, &nut, &dt, SiderealSelection::unset());
    let bound = star_bound(&star_serving, &dt, &nut);
    // Same composition as the planet serving: GAST over the
    // sealed nutation pair and the Vondrak mean obliquity.
    let jd_ut = 2460000.0;
    let jd_tt = jd_ut + 0.0008;
    let (dpsi_rad, deps_rad) = nut(jd_tt).unwrap();
    let eps_true =
        crate::positions::vondrak_mean_obliquity_deg(jd_tt).unwrap() + deps_rad.to_degrees();
    // The test Delta-T word is the constant 0.0008, so the serving's
    // construction-time offsets equal these boundary words.
    let (off_lo, off_hi) = crate::time::sidereal_continuity_offsets(0.0008, 0.0008).unwrap();
    let sdt = SiderealDeltaT::try_new(0.0008, off_lo, off_hi).unwrap();
    let expect =
        greenwich_apparent_sidereal_time_hours(jd_ut, eps_true, dpsi_rad.to_degrees(), &sdt)
            .unwrap();
    assert!((bound.gast_hours_at_ut(jd_ut).unwrap() - expect).abs() < 1e-12);
    assert!(matches!(
        bound.gast_hours_at_ut(f64::NAN),
        Err(Error::NonFiniteArgument { .. })
    ));
    assert!(matches!(
        bound.last_hours_at_ut(jd_ut, f64::INFINITY),
        Err(Error::InvalidObserver { .. })
    ));
}

#[test]
fn last_hours_matches_direct_sidtime_composition() {
    use crate::time::{SiderealDeltaT, local_apparent_sidereal_time_hours};
    let provider = Unserved;
    let nut = |_: f64| Ok(((10.0_f64 / 3_600_000.0).to_radians(), 0.0));
    let serving = SealedServing::new(&provider, &nut);
    let dt = |_: f64| Ok(0.0008);
    let bound = transit_bound(&serving, &dt);
    let jd_ut = 2460000.0;
    let lon = 12.5;
    let jd_tt = jd_ut + 0.0008;
    let (dpsi_rad, deps_rad) = nut(jd_tt).unwrap();
    let eps_true =
        crate::positions::vondrak_mean_obliquity_deg(jd_tt).unwrap() + deps_rad.to_degrees();
    let sdt = SiderealDeltaT::try_new(0.0008, bound.offset_lo_deg, bound.offset_hi_deg).unwrap();
    let direct =
        local_apparent_sidereal_time_hours(jd_ut, lon, eps_true, dpsi_rad.to_degrees(), &sdt)
            .unwrap();
    assert_eq!(bound.last_hours_at_ut(jd_ut, lon).unwrap(), direct);
}

// -- slice 4a: sealed lunar-point serving tests ----------------

fn lunar_bound<'a>(
    serving: &'a SealedServing<'a, Unserved>,
    dt: &'a dyn Fn(f64) -> Result<f64>,
    body_id: i32,
) -> SealedLunarEventServing<'a, Unserved> {
    use crate::constants::FLG_SWIEPH;
    SealedLunarEventServing::try_new(serving, dt, body_id, FLG_SWIEPH).expect("lunar id serves")
}

#[test]
fn lunar_event_body_gate_pins_served_ids_and_refusals() {
    use crate::constants::{
        CHIRON, EARTH, ECL_NUT, FLG_BARYCTR, FLG_HELCTR, FLG_J2000, FLG_NONUT, FLG_SIDEREAL,
        FLG_SWIEPH, FLG_XYZ, MEAN_APOG, MEAN_NODE, MOON, OSCU_APOG, SUN, TRUE_NODE,
    };
    use crate::constants::{INTP_APOG, INTP_PERG};
    for body in [
        MEAN_NODE, TRUE_NODE, MEAN_APOG, OSCU_APOG, INTP_APOG, INTP_PERG, -MEAN_NODE, -TRUE_NODE,
    ] {
        assert!(is_lunar_event_body(body), "body {body}");
    }
    for body in [SUN, MOON, EARTH, ECL_NUT, CHIRON, 40, 99999, -1] {
        if body == ECL_NUT {
            continue;
        }
        assert!(!is_lunar_event_body(body), "body {body}");
    }
    assert!(!is_lunar_event_body(ECL_NUT));
    let provider = Unserved;
    let nut = |_: f64| Ok((0.0, 0.0));
    let serving = SealedServing::new(&provider, &nut);
    let dt = |_: f64| Ok(0.0008);
    // Degenerate ids refuse InvalidBody.
    for body in [EARTH, ECL_NUT] {
        assert!(matches!(
            SealedLunarEventServing::try_new(&serving, &dt, body, FLG_SWIEPH),
            Err(Error::InvalidBody { body_id, .. }) if body_id == body
        ));
    }
    // Non-lunar out-of-map ids refuse UnknownBody.
    for body in [SUN, MOON, CHIRON, 17, 40, 99999] {
        assert!(matches!(
            SealedLunarEventServing::try_new(&serving, &dt, body, FLG_SWIEPH),
            Err(Error::UnknownBody { body_id }) if body_id == body
        ));
    }
    // Selections outside the verified tropical basis refuse.
    for flags in [
        FLG_SWIEPH | FLG_SIDEREAL,
        FLG_SWIEPH | FLG_J2000,
        FLG_SWIEPH | FLG_NONUT,
        FLG_SWIEPH | FLG_XYZ,
        FLG_SWIEPH | FLG_HELCTR,
        FLG_SWIEPH | FLG_BARYCTR,
    ] {
        assert!(matches!(
            SealedLunarEventServing::try_new(&serving, &dt, MEAN_NODE, flags),
            Err(Error::UnsupportedFlagCombination { .. })
        ));
    }
    // South aliases resolve to their north source.
    let south = SealedLunarEventServing::try_new(&serving, &dt, -MEAN_NODE, FLG_SWIEPH)
        .expect("south alias serves");
    assert_eq!(south.body_id(), -MEAN_NODE);
    assert_eq!(south.north_id(), MEAN_NODE);
}

#[test]
fn lunar_south_alias_is_north_antipode() {
    use crate::constants::MEAN_NODE;
    let provider = Unserved;
    let nut = |_: f64| Ok((0.0, 0.0));
    let serving = SealedServing::new(&provider, &nut);
    let dt = |_: f64| Ok(0.0008);
    // Analytic mean words need no channels, so both sides evaluate
    // over the unserved provider; the release public-entry rule is
    // the exact antipode of the north words.
    let north_bound = lunar_bound(&serving, &dt, MEAN_NODE);
    let south_bound = lunar_bound(&serving, &dt, -MEAN_NODE);
    let jd_tt = 2460000.0;
    let (nlon, nlat, ndist) = north_bound.lunar_ecliptic_at_tt(jd_tt).unwrap();
    let (slon, slat, sdist) = south_bound.lunar_ecliptic_at_tt(jd_tt).unwrap();
    assert!(((nlon + 180.0).rem_euclid(360.0) - slon).abs() < 1e-12);
    assert_eq!(slat, -nlat);
    assert_eq!(sdist, ndist);
}

#[test]
fn lunar_mean_point_search_needs_no_channels() {
    use crate::constants::{CALC_RISE, CALC_SET, MEAN_NODE};
    let provider = Unserved;
    let nut = |_: f64| Ok((0.0, 0.0));
    let serving = SealedServing::new(&provider, &nut);
    let dt = |_: f64| Ok(0.0008);
    let bound = lunar_bound(&serving, &dt, MEAN_NODE);
    // Analytic clean-room words plus the nutation stub: the search
    // runs end to end with no served channel (never a typed miss).
    let (res, tret) = lunar_rise_set_ut_serving(&bound, &rise_request(CALC_RISE)).unwrap();
    assert!(res == 0 || res == -2);
    if res == 0 {
        assert!(tret[0].is_finite() && tret[0] != 0.0);
    }
    let (res, _) = lunar_rise_set_ut_serving(&bound, &rise_request(CALC_SET)).unwrap();
    assert!(res == 0 || res == -2);
    // Transit kinds belong to the lunar transit entry.
    for kind in [
        crate::constants::CALC_MTRANSIT,
        crate::constants::CALC_ITRANSIT,
    ] {
        assert!(matches!(
            lunar_rise_set_ut_serving(&bound, &rise_request(kind)),
            Err(Error::InvalidSearchStep { .. })
        ));
    }
    // Non-finite starts refuse before any read.
    let mut bad = rise_request(CALC_RISE);
    bad.jd_start_ut = f64::NAN;
    assert!(matches!(
        lunar_rise_set_ut_serving(&bound, &bad),
        Err(Error::NonFiniteSearchValue { .. })
    ));
}

#[test]
fn lunar_true_node_propagates_sealed_miss_fail_closed() {
    use crate::constants::{CALC_RISE, TRUE_NODE};
    let provider = Unserved;
    let nut = |_: f64| Ok((0.0, 0.0));
    let serving = SealedServing::new(&provider, &nut);
    let dt = |_: f64| Ok(0.0008);
    let bound = lunar_bound(&serving, &dt, TRUE_NODE);
    // No stored channel and no Moon state: the upfront probe refuses
    // with a typed miss (the escalated support-miss naming the Moon
    // channel), never a fabricated event.
    let err = lunar_rise_set_ut_serving(&bound, &rise_request(CALC_RISE)).unwrap_err();
    assert!(matches!(
        err,
        Error::LebBodyNotFound { .. } | Error::EphemerisRange { .. } | Error::UnknownBody { .. }
    ));
    // The Hindu arm fails on the same probe, before any search.
    let mut hindu = rise_request(CALC_RISE);
    hindu.rsmi |= crate::constants::BIT_GEOCTR_NO_ECL_LAT;
    let err = lunar_rise_set_ut_serving(&bound, &hindu).unwrap_err();
    assert!(matches!(
        err,
        Error::LebBodyNotFound { .. } | Error::EphemerisRange { .. } | Error::UnknownBody { .. }
    ));
}

#[test]
fn lunar_gast_threads_nutation_and_obliquity() {
    use crate::constants::MEAN_NODE;
    use crate::time::{SiderealDeltaT, greenwich_apparent_sidereal_time_hours};
    let provider = Unserved;
    let nut = |_: f64| Ok(((10.0_f64 / 3_600_000.0).to_radians(), 0.0));
    let serving = SealedServing::new(&provider, &nut);
    let dt = |_: f64| Ok(0.0008);
    let bound = lunar_bound(&serving, &dt, MEAN_NODE);
    // Same composition as the planet serving: GAST over the
    // sealed nutation pair and the Vondrak mean obliquity.
    let jd_ut = 2460000.0;
    let jd_tt = jd_ut + 0.0008;
    let (dpsi_rad, deps_rad) = nut(jd_tt).unwrap();
    let eps_true =
        crate::positions::vondrak_mean_obliquity_deg(jd_tt).unwrap() + deps_rad.to_degrees();
    let (off_lo, off_hi) = crate::time::sidereal_continuity_offsets(0.0008, 0.0008).unwrap();
    let sdt = SiderealDeltaT::try_new(0.0008, off_lo, off_hi).unwrap();
    let expect =
        greenwich_apparent_sidereal_time_hours(jd_ut, eps_true, dpsi_rad.to_degrees(), &sdt)
            .unwrap();
    assert!((bound.gast_hours_at_ut(jd_ut).unwrap() - expect).abs() < 1e-12);
    assert!(matches!(
        bound.lunar_transit_at_ut(f64::NAN, 0.0, false),
        Err(Error::NonFiniteSearchValue { .. })
    ));
}

// -- slice 4b: sealed minor-body serving tests ----------------

fn minor_bound<'a>(
    serving: &'a SealedServing<'a, Unserved>,
    dt: &'a dyn Fn(f64) -> Result<f64>,
    body_id: i32,
) -> SealedMinorEventServing<'a, Unserved> {
    use crate::constants::FLG_SWIEPH;
    SealedMinorEventServing::try_new(serving, dt, body_id, FLG_SWIEPH).expect("minor id serves")
}

#[test]
fn minor_event_body_gate_pins_served_ids_and_refusals() {
    use crate::constants::{
        CERES, CHIRON, EARTH, ECL_NUT, FLG_BARYCTR, FLG_HELCTR, FLG_J2000, FLG_NONUT, FLG_SIDEREAL,
        FLG_SWIEPH, FLG_XYZ, JUNO, MEAN_NODE, MOON, PALLAS, PHOLUS, SUN, VESTA,
    };
    // Dedicated ids, AST_OFFSET aliases and direct numbers serve.
    for body in [
        CHIRON, PHOLUS, CERES, PALLAS, JUNO, VESTA, 10001, 10002, 10003, 10004, 12060, 15145,
        10433, 11181, 111955, 146199, 172173, 235088,
    ] {
        assert!(is_minor_event_body(body), "body {body}");
    }
    // Planets, lunar points, fictitious bodies, moons, bogus ids and
    // non-table numbers refuse (their bindings own them, or the
    // release probe maps them to the illegal-body refusal).
    for body in [
        SUN, MOON, 3, MEAN_NODE, 40, 99999, 109999, 10015, 1000001, -10, -11,
    ] {
        assert!(!is_minor_event_body(body), "body {body}");
    }
    assert!(!is_minor_event_body(EARTH));
    assert!(!is_minor_event_body(ECL_NUT));
    let provider = Unserved;
    let nut = |_: f64| Ok((0.0, 0.0));
    let serving = SealedServing::new(&provider, &nut);
    let dt = |_: f64| Ok(0.0008);
    // Degenerate ids refuse InvalidBody.
    for body in [EARTH, ECL_NUT] {
        assert!(matches!(
            SealedMinorEventServing::try_new(&serving, &dt, body, FLG_SWIEPH),
            Err(Error::InvalidBody { body_id, .. }) if body_id == body
        ));
    }
    // Non-table ids refuse UnknownBody.
    for body in [SUN, MOON, MEAN_NODE, 40, 99999, 109999] {
        assert!(matches!(
            SealedMinorEventServing::try_new(&serving, &dt, body, FLG_SWIEPH),
            Err(Error::UnknownBody { body_id }) if body_id == body
        ));
    }
    // Selections outside the verified tropical basis refuse.
    for flags in [
        FLG_SWIEPH | FLG_SIDEREAL,
        FLG_SWIEPH | FLG_J2000,
        FLG_SWIEPH | FLG_NONUT,
        FLG_SWIEPH | FLG_XYZ,
        FLG_SWIEPH | FLG_HELCTR,
        FLG_SWIEPH | FLG_BARYCTR,
    ] {
        assert!(matches!(
            SealedMinorEventServing::try_new(&serving, &dt, CHIRON, flags),
            Err(Error::UnsupportedFlagCombination { .. })
        ));
    }
    // Aliases resolve to their dedicated table id; the unserved stub
    // channels nothing, so every table id takes the fallback arm.
    let aliased = minor_bound(&serving, &dt, 10001);
    assert_eq!(aliased.body_id(), 10001);
    assert_eq!(aliased.minor_id(), CERES);
    assert!(!aliased.is_leb_channeled());
    let direct = minor_bound(&serving, &dt, CHIRON);
    assert_eq!(direct.minor_id(), CHIRON);
}

#[test]
fn minor_table_key_space_matches_release_inventory() {
    // The gate is table membership, so the catalogue table must carry
    // exactly the 37 release `MINOR_BODY_ELEMENTS` keys (dedicated
    // 15-20 plus the 31 `AST_OFFSET + N` numbers, from the
    // pinned release): any drift silently re-gates the event path.
    let mut ported: Vec<i32> = crate::minor_body_tables::MINOR_BODY_TABLE
        .iter()
        .map(|row| row.body_id)
        .collect();
    ported.sort_unstable();
    let release_keys = [
        15, 16, 17, 18, 19, 20, 10010, 10016, 10052, 10055, 10080, 10087, 10433, 10511, 10704,
        10944, 11181, 11221, 11566, 11685, 14179, 17066, 18405, 20199, 30000, 35143, 38978, 60000,
        100377, 100482, 109942, 111955, 146108, 146199, 146472, 172173, 235088,
    ];
    assert_eq!(ported, release_keys);
    for body in release_keys {
        assert!(is_minor_event_body(body), "key {body}");
    }
}

#[test]
fn minor_fallback_propagates_sealed_miss_fail_closed() {
    use crate::constants::{CALC_RISE, FLG_SWIEPH};
    // Bennu (111955) is the table id without a base LEB channel, so
    // the unserved stub routes it to the fallback composition, whose
    // Earth read fails with the typed miss (escalated), never a
    // fabricated event.
    let provider = Unserved;
    let nut = |_: f64| Ok((0.0, 0.0));
    let serving = SealedServing::new(&provider, &nut);
    let dt = |_: f64| Ok(0.0008);
    let bound =
        SealedMinorEventServing::try_new(&serving, &dt, 111955, FLG_SWIEPH).expect("bennu gates");
    assert!(!bound.is_leb_channeled());
    let err = minor_rise_set_ut_serving(&bound, &rise_request(CALC_RISE)).unwrap_err();
    assert!(matches!(
        err,
        Error::LebBodyNotFound { .. } | Error::EphemerisRange { .. } | Error::UnknownBody { .. }
    ));
    // The transit arm fails on the same probe, before any search.
    let err = minor_transit_ut_serving(&bound, 2460400.0, 12.5, false).unwrap_err();
    assert!(matches!(
        err,
        Error::LebBodyNotFound { .. } | Error::EphemerisRange { .. } | Error::UnknownBody { .. }
    ));
    // Non-finite epochs refuse before any read.
    assert!(matches!(
        bound.kepler_ecliptic_at_tt(f64::NAN),
        Err(Error::NonFiniteArgument { .. })
    ));
    assert!(matches!(
        bound.kepler_transit_at_ut(f64::NAN, 0.0, false),
        Err(Error::NonFiniteSearchValue { .. })
    ));
}

#[test]
fn minor_gast_threads_nutation_and_obliquity() {
    use crate::constants::CHIRON;
    use crate::time::{SiderealDeltaT, greenwich_apparent_sidereal_time_hours};
    let provider = Unserved;
    let nut = |_: f64| Ok(((10.0_f64 / 3_600_000.0).to_radians(), 0.0));
    let serving = SealedServing::new(&provider, &nut);
    let dt = |_: f64| Ok(0.0008);
    let bound = minor_bound(&serving, &dt, CHIRON);
    // Same composition as the planet serving: GAST over the
    // sealed nutation pair and the Vondrak mean obliquity.
    let jd_ut = 2460000.0;
    let jd_tt = jd_ut + 0.0008;
    let (dpsi_rad, deps_rad) = nut(jd_tt).unwrap();
    let eps_true =
        crate::positions::vondrak_mean_obliquity_deg(jd_tt).unwrap() + deps_rad.to_degrees();
    let (off_lo, off_hi) = crate::time::sidereal_continuity_offsets(0.0008, 0.0008).unwrap();
    let sdt = SiderealDeltaT::try_new(0.0008, off_lo, off_hi).unwrap();
    let expect =
        greenwich_apparent_sidereal_time_hours(jd_ut, eps_true, dpsi_rad.to_degrees(), &sdt)
            .unwrap();
    assert!((bound.gast_hours_at_ut(jd_ut).unwrap() - expect).abs() < 1e-12);
}

// -- slice 4c: sealed fictitious-body serving tests --------------

fn fict_bound<'a>(
    serving: &'a SealedServing<'a, Unserved>,
    dt: &'a dyn Fn(f64) -> Result<f64>,
    body_id: i32,
) -> SealedFictitiousEventServing<'a, Unserved> {
    use crate::constants::FLG_SWIEPH;
    SealedFictitiousEventServing::try_new(serving, dt, body_id, FLG_SWIEPH)
        .expect("fictitious id serves")
}

#[test]
fn fict_event_body_gate_pins_served_ids_and_refusals() {
    use crate::constants::{
        CHIRON, CUPIDO, EARTH, ECL_NUT, FLG_BARYCTR, FLG_HELCTR, FLG_J2000, FLG_NONUT, FLG_RADIANS,
        FLG_SIDEREAL, FLG_SWIEPH, FLG_XYZ, HARRINGTON, ISIS, MEAN_NODE, MOON, NEPTUNE_ADAMS,
        NEPTUNE_LEVERRIER, NIBIRU, PLUTO_LOWELL, PLUTO_PICKERING, POSEIDON, PROSERPINA, SUN,
        VULCAN, WALDEMATH, WHITE_MOON,
    };
    // All 18 model-served ids gate open (Uranian 40-47, Isis 48,
    // predictions/Vulcan/Proserpina 50-55/57, White Moon 56,
    // Waldemath 58).
    for body in [
        CUPIDO,
        41,
        42,
        43,
        44,
        45,
        46,
        POSEIDON,
        ISIS,
        HARRINGTON,
        NEPTUNE_LEVERRIER,
        NEPTUNE_ADAMS,
        PLUTO_LOWELL,
        PLUTO_PICKERING,
        VULCAN,
        WHITE_MOON,
        PROSERPINA,
        WALDEMATH,
    ] {
        assert!(is_fictitious_event_body(body), "body {body}");
    }
    // Planets, lunar points, minor bodies, model-less fictitious ids
    // and bogus ids refuse (their bindings own them, or the release
    // probe maps them to the illegal-body refusal).
    for body in [
        SUN, MOON, 3, MEAN_NODE, CHIRON, 10001, NIBIRU, 59, 69, 70, 99999, -1,
    ] {
        assert!(!is_fictitious_event_body(body), "body {body}");
    }
    assert!(!is_fictitious_event_body(EARTH));
    assert!(!is_fictitious_event_body(ECL_NUT));
    let provider = Unserved;
    let nut = |_: f64| Ok((0.0, 0.0));
    let serving = SealedServing::new(&provider, &nut);
    let dt = |_: f64| Ok(0.0008);
    // Degenerate ids refuse InvalidBody.
    for body in [EARTH, ECL_NUT] {
        assert!(matches!(
            SealedFictitiousEventServing::try_new(&serving, &dt, body, FLG_SWIEPH),
            Err(Error::InvalidBody { body_id, .. }) if body_id == body
        ));
    }
    // Model-less recognised ids and non-hypothetical ids refuse
    // UnknownBody.
    for body in [NIBIRU, 59, 69, SUN, MOON, MEAN_NODE, CHIRON, 99999] {
        assert!(matches!(
            SealedFictitiousEventServing::try_new(&serving, &dt, body, FLG_SWIEPH),
            Err(Error::UnknownBody { body_id }) if body_id == body
        ));
    }
    // Selections outside the verified tropical basis refuse.
    for flags in [
        FLG_SWIEPH | FLG_SIDEREAL,
        FLG_SWIEPH | FLG_J2000,
        FLG_SWIEPH | FLG_NONUT,
        FLG_SWIEPH | FLG_XYZ,
        FLG_SWIEPH | FLG_RADIANS,
        FLG_SWIEPH | FLG_HELCTR,
        FLG_SWIEPH | FLG_BARYCTR,
    ] {
        assert!(matches!(
            SealedFictitiousEventServing::try_new(&serving, &dt, CUPIDO, flags),
            Err(Error::UnsupportedFlagCombination { .. })
        ));
    }
    // Served construction carries the requested id and flags.
    let bound = fict_bound(&serving, &dt, WHITE_MOON);
    assert_eq!(bound.body_id(), WHITE_MOON);
    assert_eq!(bound.calc_flags(), FLG_SWIEPH);
}

#[test]
// Keep the recorded source fixture literals and their exact binary64 words.
#[allow(clippy::excessive_precision)]
fn fict_aberration_helpers_match_independent_pins() {
    // Independent pins evaluated from the release formulas through a
    // fresh interpreter session (slice-4c work notes); bitwise pins
    // catch association/constant slips in the port.
    let geo = [1.0, 0.5, 0.25];
    let vel = [0.0, 0.0172, 0.001];
    let b = bradley_aberration_xyz(geo, vel);
    assert_eq!(b[0].to_bits(), 0.99995538066620093_f64.to_bits());
    assert_eq!(b[1].to_bits(), 0.50009149735734859_f64.to_bits());
    assert_eq!(b[2].to_bits(), 0.24999546185400651_f64.to_bits());
    let rg = (geo[0] * geo[0] + geo[1] * geo[1] + geo[2] * geo[2]).sqrt();
    let f = first_order_aberration_xyz(geo, vel, rg);
    assert_eq!(f[0].to_bits(), 1.0_f64.to_bits());
    assert_eq!(f[1].to_bits(), 0.50011380702469821_f64.to_bits());
    assert_eq!(f[2].to_bits(), 0.25000661668748247_f64.to_bits());
    let s = ecl_xyz_to_spherical([0.3, -1.2, 4.5]);
    assert_eq!(s.0.to_bits(), 284.03624346792645_f64.to_bits());
    assert_eq!(s.1.to_bits(), 74.630476198695646_f64.to_bits());
    assert_eq!(s.2.to_bits(), 4.6669047558312142_f64.to_bits());
    // Degenerate shaping arm: the zero vector reports zero latitude.
    assert_eq!(ecl_xyz_to_spherical([0.0, 0.0, 0.0]), (0.0, 0.0, 0.0));
}

#[test]
fn fict_bradley_aberration_properties() {
    // Zero observer velocity is the identity, bitwise.
    let geo = [1.0, 0.5, 0.25];
    assert_eq!(bradley_aberration_xyz(geo, [0.0, 0.0, 0.0]), geo);
    // Sub-floor inputs pass through unchanged.
    assert_eq!(
        bradley_aberration_xyz([1e-31, 0.0, 0.0], [0.01, 0.0, 0.0]),
        [1e-31, 0.0, 0.0]
    );
    // Renormalization preserves the distance (the corrected direction
    // is re-scaled to the input range).
    let out = bradley_aberration_xyz(geo, [0.0, 0.0172, 0.001]);
    let want = (geo[0] * geo[0] + geo[1] * geo[1] + geo[2] * geo[2]).sqrt();
    let got = (out[0] * out[0] + out[1] * out[1] + out[2] * out[2]).sqrt();
    assert!((got - want).abs() < 1e-15 * want);
    // The shift points along the observer velocity.
    let dot = (out[0] - geo[0]) * 0.0 + (out[1] - geo[1]) * 0.0172 + (out[2] - geo[2]) * 0.001;
    assert!(dot > 0.0);
    // First-order shift magnitude is |v| * rg / c.
    let rg = want;
    let f = first_order_aberration_xyz(geo, [0.0, 0.0172, 0.001], rg);
    let shift =
        ((f[0] - geo[0]).powi(2) + (f[1] - geo[1]).powi(2) + (f[2] - geo[2]).powi(2)).sqrt();
    let vmag = (0.0172_f64.powi(2) + 0.001_f64.powi(2)).sqrt();
    // Norm rounding after the per-component shift: a few ulps of the
    // ~1e-4 shift (measured 2e-17 on this machine).
    assert!((shift - vmag * rg / crate::positions::C_LIGHT_AU_DAY).abs() < 1e-16);
}

#[test]
fn fict_rotate_to_of_date_preserves_magnitude_and_j2000_identity() {
    // At J2000 the precession early-out leaves the direction
    // untouched; only the lon/lat→xyz round-trip moves words (a few
    // ulps — measured 4e-19 here).
    let v = [0.001, -0.002, 0.0005];
    let w0 = rotate_vector_to_of_date(v, crate::hypothetical::J2000_JD).unwrap();
    for i in 0..3 {
        assert!((w0[i] - v[i]).abs() < 1e-15, "lane {i}");
    }
    // Elsewhere the magnitude is preserved (direction precessed).
    let w = rotate_vector_to_of_date(v, 2460000.0).unwrap();
    let mag_v = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    let mag_w = (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt();
    assert!((mag_w - mag_v).abs() < 1e-18);
    assert!(w != v);
    // The zero vector maps to zeros.
    assert_eq!(
        rotate_vector_to_of_date([0.0, 0.0, 0.0], 2460000.0).unwrap(),
        [0.0, 0.0, 0.0]
    );
}

#[test]
fn fict_white_moon_model_words_are_linear_and_geocentric() {
    // BODY-03 symbolic words: uniform longitude rate, zero latitude,
    // positive distance (composition detail stays sweep-verified).
    let m0 = crate::hypothetical::calc_white_moon_position(2451545.0, false).unwrap();
    let m1 = crate::hypothetical::calc_white_moon_position(2451645.0, false).unwrap();
    assert_eq!(m0.lat_deg, 0.0);
    assert!(m0.dist_au > 0.0);
    let rate = crate::hypothetical::white_moon_rate_deg_per_day();
    assert!(((m1.lon_deg - m0.lon_deg).rem_euclid(360.0) - rate * 100.0).abs() < 1e-9);
}

#[test]
fn fict_serving_propagates_sealed_miss_fail_closed() {
    use crate::constants::{CALC_RISE, CUPIDO, FLG_SWIEPH};
    // Cupido on the unserved stub: the upfront probe fails on the
    // Earth/Sun read with the typed miss (escalated), never a
    // fabricated event.
    let provider = Unserved;
    let nut = |_: f64| Ok((0.0, 0.0));
    let serving = SealedServing::new(&provider, &nut);
    let dt = |_: f64| Ok(0.0008);
    let bound = SealedFictitiousEventServing::try_new(&serving, &dt, CUPIDO, FLG_SWIEPH)
        .expect("cupido gates");
    let err = fictitious_rise_set_ut_serving(&bound, &rise_request(CALC_RISE)).unwrap_err();
    assert!(matches!(
        err,
        Error::LebBodyNotFound { .. } | Error::EphemerisRange { .. } | Error::UnknownBody { .. }
    ));
    // The transit arm fails on the same probe, before any search.
    let err = fictitious_transit_ut_serving(&bound, 2460400.0, 12.5, false).unwrap_err();
    assert!(matches!(
        err,
        Error::LebBodyNotFound { .. } | Error::EphemerisRange { .. } | Error::UnknownBody { .. }
    ));
    // Non-finite epochs refuse before any read.
    assert!(matches!(
        bound.fict_geo_ecliptic_at_tt(f64::NAN),
        Err(Error::NonFiniteArgument { .. })
    ));
    assert!(matches!(
        bound.fict_transit_at_ut(f64::NAN, 0.0, false),
        Err(Error::NonFiniteSearchValue { .. })
    ));
    assert!(matches!(
        bound.fict_topo_ecliptic_at_ut(f64::NAN, 0.0, 0.0, 0.0),
        Err(Error::NonFiniteArgument { .. })
    ));
}
