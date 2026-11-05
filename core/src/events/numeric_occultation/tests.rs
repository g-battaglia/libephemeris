//! Pure callback ordering/ownership controls. Primary compiles and runs these.
use super::*;
use crate::events::{ECL_RSUN_AU, eclipse_where_core_via, sol_how_core_via};

fn geometry(quiet: bool) -> WhereGeometry {
    eclipse_where_core_via(
        if quiet {
            (0.0, 0.0025, 0.0)
        } else {
            (0.0025, 0.0, 0.0)
        },
        (1.0, 0.0, 0.0),
        17.0,
        ECL_RSUN_AU,
    )
    .unwrap()
}

fn circumstances() -> SolHowCircumstances {
    sol_how_core_via(
        (0.0, 0.0, 1.0),
        (0.0, 0.0, 0.0025),
        ECL_RSUN_AU,
        true,
        true,
        (12.0, -30.0, -29.0),
        0.0,
        (1.0, 2.0),
    )
    .unwrap()
}

struct Source {
    quiet: bool,
    global_error: Option<Error>,
    local_error: Option<Error>,
    seen: Vec<(&'static str, u64, i32, i32, [f64; 3])>,
}
impl Source {
    fn new(quiet: bool) -> Self {
        Self {
            quiet,
            global_error: None,
            local_error: None,
            seen: Vec::new(),
        }
    }
}
impl NumericOccultWhereSource for Source {
    fn global(&mut self, jd: f64, body: i32, flags: i32) -> Result<WhereGeometry, Error> {
        self.seen
            .push(("global", jd.to_bits(), body, flags, [0.0; 3]));
        if let Some(e) = &self.global_error {
            return Err(e.clone());
        }
        Ok(geometry(self.quiet))
    }
    fn local(
        &mut self,
        jd: f64,
        body: i32,
        flags: i32,
        site: [f64; 3],
    ) -> Result<SolHowCircumstances, Error> {
        self.seen.push(("local", jd.to_bits(), body, flags, site));
        if let Some(e) = &self.local_error {
            return Err(e.clone());
        }
        Ok(circumstances())
    }
}

#[test]
fn geometry_then_local_keep_epoch_target_flags_and_actual_center() {
    let mut source = Source::new(false);
    let jd = 2451545.0;
    let flags = 0x8002;
    let result = numeric_occult_where_via(&mut source, jd, 3, flags).unwrap();
    let g = geometry(false);
    assert_eq!(
        source.seen,
        vec![
            ("global", jd.to_bits(), 3, flags, [0.0; 3]),
            (
                "local",
                jd.to_bits(),
                3,
                flags,
                [g.lon_deg(), g.lat_deg(), 0.0]
            ),
        ]
    );
    assert_eq!(result.retflag(), g.retc());
    assert_eq!(result.geopos()[2..], [0.0; 8]);
    assert_eq!(result.attr()[3].to_bits(), g.dcore()[0].to_bits());
    assert_eq!(result.attr()[6], -29.0);
}

#[test]
fn zero_global_character_still_consumes_local_and_keeps_closest_point() {
    let mut source = Source::new(true);
    let result = numeric_occult_where_via(&mut source, 2451545.0, 3, 2).unwrap();
    assert_eq!(geometry(true).retc(), 0);
    assert_eq!(result.retflag(), 0);
    assert_eq!(source.seen.len(), 2);
    assert_eq!(result.geopos()[0], geometry(true).lon_deg());
    assert_eq!(result.attr()[4], 12.0);
}

#[test]
fn global_error_prevents_local_read_and_preserves_owned_carrier() {
    let mut source = Source::new(false);
    let marker = Error::UnknownBody { body_id: 333 };
    source.global_error = Some(marker.clone());
    assert_eq!(
        numeric_occult_where_via(&mut source, 2451545.0, 3, 2),
        Err(marker)
    );
    assert_eq!(source.seen.len(), 1);
    assert_eq!(source.seen[0].1, 2451545.0_f64.to_bits());
}

#[test]
fn local_error_is_not_hidden_by_global_zero_character() {
    let mut source = Source::new(true);
    let marker = Error::UnknownBody { body_id: 334 };
    source.local_error = Some(marker.clone());
    assert_eq!(
        numeric_occult_where_via(&mut source, 1.0, 3, 2),
        Err(marker)
    );
    assert_eq!(source.seen.len(), 2);
}

#[test]
fn target_radius_domain_is_checked_before_unneeded_horizon_geometry() {
    let error = numeric_occult_how_via(
        (0.0, 0.0, 1.0 / 1.1),
        (0.0, 0.0, 0.0),
        1.0,
        false,
        &mut || panic!("unused horizon"),
        0.0,
        &mut || panic!("unused Saros"),
    )
    .unwrap_err();
    assert!(matches!(error, Error::MathDomain { value } if (value - 1.1).abs() < 1e-15));
}

#[test]
fn zero_distance_and_impossible_lunar_radius_have_typed_errors() {
    for radius in [0.0, 1.0] {
        let error = numeric_occult_how_via(
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 0.0),
            radius,
            false,
            &mut || panic!("unused horizon"),
            0.0,
            &mut || panic!("unused Saros"),
        )
        .unwrap_err();
        assert_eq!(error, Error::RationalDivisionByZero);
    }
    let error = numeric_occult_how_via(
        (0.0, 0.0, 1.0),
        (0.0, 0.0, super::super::ECL_RMOON_AU / 2.0),
        0.0,
        false,
        &mut || panic!("unused horizon"),
        0.0,
        &mut || panic!("unused Saros"),
    )
    .unwrap_err();
    assert!(matches!(error, Error::MathDomain { value } if value == 2.0));
}

#[test]
fn coincident_discs_cover_the_target_and_keep_the_circumstance_layout() {
    let value = numeric_occult_how_via(
        (0.0, 0.0, 1.0),
        (0.0, 0.0, 0.0023),
        ECL_RSUN_AU,
        true,
        &mut || Ok((180.0, 60.0, 60.01)),
        0.0,
        &mut || Ok((1.0, 2.0)),
    )
    .unwrap();
    let attributes = value.attr();
    assert!(attributes[0] > 0.0);
    assert!(attributes[2] >= 1.0);
    assert_eq!(attributes[4], 180.0);
    assert_eq!(attributes[5], 60.0);
    assert_eq!(attributes[9], 1.0);
    assert_eq!(attributes[10], 2.0);
}
