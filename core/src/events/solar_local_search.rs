//! Native public local solar search policy over common contact kernels.
//! Provider, horizon and circumstance calls retain their ordered errors; only
//! the global-candidate catch has a supplied, constructor-aware predicate.
use super::solar_eclipses::bisect_sign_change;
use super::*;

/// Search at original UT1 JD for a visible local eclipse. `eph` owns its site
/// and acquisition channels. `backwards` selects the candidate-collection arm.
/// Global candidate errors satisfying `caught_global` break collection or raise
/// the provided `exhausted` constructor; other failures are terminal. Contact,
/// maximum, visibility and horizon errors propagate with no range conversion.
/// Returns character bits, ten UT1 JDs and twenty circumstances; absent contacts
/// and the final three time slots are zero. Uses 120 candidate iterations, a
/// golden-section maximum resolved to the input clock's representable spacing,
/// and existing 60-step sign-change contact bisections.
pub fn sol_eclipse_when_loc_source_via(
    eph: &impl SolarLocalEphemeris,
    tjdut: f64,
    backwards: bool,
    caught_global: &dyn Fn(&Error) -> bool,
    exhausted: &dyn Fn(Option<Error>) -> Error,
) -> Result<(i32, [f64; 10], [f64; 20])> {
    let mut jd = tjdut;
    for _ in 0..SOL_LOC_MAX_ECLIPSES {
        // Forward candidates start0.35day before every iteration; local
        // contact/maximum epoch gates decide whether that candidate is kept.
        // Backward: collect forward from a year back and keep the
        // most recent eclipse before the epoch window.
        let jd_max_global = if backwards {
            let earlier_jd = jd - SOL_LOC_BACKWARD_COLLECT_DAYS;
            let horizon = jd + WHEN_LOC_EPOCH_WINDOW_DAYS;
            let mut found: Option<(i32, [f64; 10])> = None;
            let mut temp_jd = earlier_jd;
            while temp_jd < horizon {
                match eph.next_global_forward(temp_jd) {
                    Ok((global_type, global_times)) => {
                        if global_times[0] < horizon {
                            found = Some((global_type, global_times));
                        }
                        temp_jd = global_times[0] + ECLIPSE_WALK_ADVANCE_DAYS;
                    }
                    Err(error) if caught_global(&error) => break,
                    Err(error) => return Err(error),
                }
            }
            match found {
                Some((_, global_times)) => global_times[0],
                None => {
                    jd -= SOL_LOC_BACKWARD_COLLECT_DAYS;
                    continue;
                }
            }
        } else {
            let anchor = jd - WHEN_LOC_EPOCH_WINDOW_DAYS;
            match eph.next_global_forward(anchor) {
                Ok((_, global_times)) => global_times[0],
                Err(error) if caught_global(&error) => return Err(exhausted(Some(error))),
                Err(error) => return Err(error),
            }
        };

        let mut separation = |jd: f64| eph.separation_deg(jd);
        let jd_local_max = find_local_maximum_source_via(
            &mut separation,
            jd_max_global,
            LOCAL_MAX_SEARCH_RANGE_DAYS,
        )?;
        let (rsun, rmoon, min_separation) = eph.radii_sep_deg(jd_local_max, false)?;
        if min_separation > rsun + rmoon {
            jd = if backwards {
                jd_max_global - 1.0
            } else {
                jd_max_global + ECLIPSE_WALK_ADVANCE_DAYS
            };
            continue;
        }
        // Local phase at maximum (ordering: the
        // `rsun - rmoon` arm can only fire for an annular-geometry
        // Moon; a negative difference falls through to the `abs` arm).
        let phase = if min_separation < rsun - rmoon {
            ECL_ANNULAR
        } else if min_separation < (rsun - rmoon).abs() {
            ECL_TOTAL
        } else {
            ECL_PARTIAL
        };
        let two_hours = 2.0 / 24.0;
        let mut f_outer = |jd_c: f64| -> Result<f64> {
            let (rs_c, rm_c, sep_c) = eph.radii_sep_deg(jd_c, false)?;
            Ok(sep_c - (rs_c + rm_c))
        };
        let mut f_inner = |jd_c: f64| -> Result<f64> {
            let (rs_c, rm_c, sep_c) = eph.radii_sep_deg(jd_c, true)?;
            Ok((rs_c - rm_c).abs() - sep_c)
        };
        let mut jd_first =
            bisect_sign_change(&mut f_outer, jd_local_max - two_hours, jd_local_max)?;
        if jd_first == 0.0 {
            jd_first =
                bisect_sign_change(&mut f_outer, jd_local_max - 2.0 * two_hours, jd_local_max)?;
        }
        let mut jd_fourth =
            bisect_sign_change(&mut f_outer, jd_local_max, jd_local_max + two_hours)?;
        if jd_fourth == 0.0 {
            jd_fourth =
                bisect_sign_change(&mut f_outer, jd_local_max, jd_local_max + 2.0 * two_hours)?;
        }
        // Epoch gate on the LOCAL contacts: an
        // eclipse whose local phases are still unfolding at the
        // epoch is kept; a wholly past (forward) or wholly future
        // (backward) one is skipped. The second arm stops the
        // re-entrant `jd = tret[0]; when_loc(jd)` idiom from
        // re-finding the same eclipse.
        let first_contact = if jd_first != 0.0 {
            jd_first
        } else {
            jd_local_max
        };
        let last_contact = if jd_fourth != 0.0 {
            jd_fourth
        } else {
            jd_local_max
        };
        if (!backwards
            && (last_contact <= tjdut + ECLIPSE_WHEN_EPOCH_MARGIN_DAYS
                || jd_local_max <= tjdut + ECLIPSE_WHEN_EPOCH_MARGIN_DAYS))
            || (backwards
                && (first_contact >= tjdut - ECLIPSE_WHEN_EPOCH_MARGIN_DAYS
                    || jd_local_max >= tjdut - ECLIPSE_WHEN_EPOCH_MARGIN_DAYS))
        {
            jd = if backwards {
                jd_max_global - 1.0
            } else {
                jd_max_global + ECLIPSE_WALK_ADVANCE_DAYS
            };
            continue;
        }
        let mut jd_second = 0.0;
        let mut jd_third = 0.0;
        if min_separation < (rsun - rmoon).abs() {
            let half_hour = LOCAL_INNER_CONTACT_SEARCH_DAYS;
            jd_second = bisect_sign_change(&mut f_inner, jd_local_max - half_hour, jd_local_max)?;
            jd_third = bisect_sign_change(&mut f_inner, jd_local_max, jd_local_max + half_hour)?;
        }
        // Visibility bits: a phase counts when the
        // Sun's APPARENT altitude is positive there; `tret` keeps the
        // geometric times regardless.
        let mut ecl_type = phase;
        for (tc, bit) in [
            (jd_local_max, ECL_MAX_VISIBLE),
            (jd_first, ECL_1ST_VISIBLE),
            (jd_second, ECL_2ND_VISIBLE),
            (jd_third, ECL_3RD_VISIBLE),
            (jd_fourth, ECL_4TH_VISIBLE),
        ] {
            if tc == 0.0 {
                continue;
            }
            if eph.sun_altaz_deg(tc)?.2 > 0.0 {
                ecl_type |= ECL_VISIBLE | bit;
            }
        }
        if ecl_type & ECL_VISIBLE == 0 {
            jd = if backwards {
                jd_max_global - 1.0
            } else {
                jd_max_global + ECLIPSE_WALK_ADVANCE_DAYS
            };
            continue;
        }
        // Sunrise/sunset bounding the visible part: a circumpolar Sun (`-2`) leaves the geometric
        // times untouched; a wholly nocturnal eclipse advances the
        // walk; a horizon crossing inside the contacts is reported
        // and, when the geometric maximum itself is down, becomes
        // the reported maximum.
        let mut jd_reported_max = jd_local_max;
        let mut jd_sunrise = 0.0;
        let mut jd_sunset = 0.0;
        let mut clamped = false;
        let crossings = eph.sun_horizon_crossings(first_contact - 0.001)?;
        let circumpolar = crossings.rise_rc == -2 || crossings.set_rc == -2;
        if !circumpolar {
            let tjdr = crossings.rise_jd;
            let tjds = crossings.set_jd;
            if tjds < jd_first || (tjds > tjdr && tjdr > jd_fourth) {
                jd = if backwards {
                    jd_max_global - 1.0
                } else {
                    jd_max_global + ECLIPSE_WALK_ADVANCE_DAYS
                };
                continue;
            }
            if jd_first < tjdr && tjdr < jd_fourth {
                jd_sunrise = tjdr;
                if ecl_type & ECL_MAX_VISIBLE == 0 {
                    jd_reported_max = tjdr;
                    clamped = true;
                }
            }
            if jd_first < tjds && tjds < jd_fourth {
                jd_sunset = tjds;
                if ecl_type & ECL_MAX_VISIBLE == 0 {
                    jd_reported_max = tjds;
                    clamped = true;
                }
            }
        }
        let (retc_how, mut attr) = eph.how_core(jd_reported_max)?;
        if clamped {
            ecl_type &= !(ECL_TOTAL | ECL_ANNULAR | ECL_PARTIAL);
            ecl_type |= retc_how & (ECL_TOTAL | ECL_ANNULAR | ECL_PARTIAL);
        }
        let (rc_where, dcore) = eph.where_core(jd_reported_max)?;
        ecl_type |= rc_where & ECL_NONCENTRAL;
        attr[3] = dcore[0];
        let tret = [
            jd_reported_max,
            jd_first,
            jd_second,
            jd_third,
            jd_fourth,
            jd_sunrise,
            jd_sunset,
            0.0,
            0.0,
            0.0,
        ];
        return Ok((ecl_type, tret, attr));
    }
    Err(exhausted(None))
}

/// Separation minimum in a symmetric UT1-day interval.
///
/// Resolve the bracket down to distinct representable epochs: a fixed
/// 1e-7-day stop can leave millisecond timing errors that become significant
/// angular circumstance errors near the zenith. At most eighty contractions
/// cover the public three-hour half-window at ordinary Julian-day precision.
/// Provider errors propagate; clock resolution is not a physical error bound.
fn find_local_maximum_source_via(
    separation: &mut impl FnMut(f64) -> Result<f64>,
    center: f64,
    range: f64,
) -> Result<f64> {
    super::solar_eclipses::find_local_maximum_via(separation, center, range)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::RefCell;

    #[test]
    fn local_maximum_resolves_submillisecond_circumstances() {
        let center = 2_451_545.0;
        let minimum = center + 1.0 / 73.0;
        let mut separation = |jd: f64| Ok(0.001_f64.hypot(12.0 * (jd - minimum)));
        let result = find_local_maximum_source_via(&mut separation, center, 0.125).unwrap();
        assert!((result - minimum).abs() <= minimum.next_up() - minimum);
    }

    struct Toy {
        calls: RefCell<Vec<f64>>,
        global: Error,
        local: bool,
    }
    impl SolarLocalEphemeris for Toy {
        fn next_global_forward(&self, jd: f64) -> Result<(i32, [f64; 10])> {
            self.calls.borrow_mut().push(jd);
            if self.local {
                let mut t = [0.0; 10];
                t[0] = 10.0;
                Ok((0, t))
            } else {
                Err(self.global.clone())
            }
        }
        fn separation_deg(&self, _: f64) -> Result<f64> {
            Err(self.global.clone())
        }
        fn radii_sep_deg(&self, _: f64, _: bool) -> Result<(f64, f64, f64)> {
            unreachable!()
        }
        fn sun_altaz_deg(&self, _: f64) -> Result<(f64, f64, f64)> {
            unreachable!()
        }
        fn topo_diameters_deg(&self, _: f64) -> Result<(f64, f64)> {
            unreachable!()
        }
        fn sun_horizon_crossings(&self, _: f64) -> Result<BodyHorizonCrossings> {
            unreachable!()
        }
        fn how_core(&self, _: f64) -> Result<(i32, [f64; 20])> {
            unreachable!()
        }
        fn where_core(&self, _: f64) -> Result<(i32, [f64; 10])> {
            unreachable!()
        }
    }
    #[test]
    fn raw_epoch_reaches_global_callback_before_finite_admission() {
        let toy = Toy {
            calls: RefCell::new(Vec::new()),
            global: Error::UnknownBody { body_id: 777 },
            local: false,
        };
        let result = sol_eclipse_when_loc_source_via(&toy, f64::NAN, false, &|_| false, &|_| {
            panic!("unrelated global failure must propagate")
        });
        assert_eq!(result, Err(Error::UnknownBody { body_id: 777 }));
        assert_eq!(toy.calls.borrow().len(), 1);
        assert!(toy.calls.borrow()[0].is_nan());
    }
    #[test]
    fn caught_forward_failure_has_active_error_but_local_range_does_not() {
        let marker = Error::UnknownBody { body_id: 777 };
        let toy = Toy {
            calls: RefCell::new(Vec::new()),
            global: marker.clone(),
            local: false,
        };
        let active = RefCell::new(None);
        let result = sol_eclipse_when_loc_source_via(&toy, 9.0, false, &|_| true, &|error| {
            *active.borrow_mut() = error;
            Error::UnknownBody { body_id: 778 }
        });
        assert_eq!(result, Err(Error::UnknownBody { body_id: 778 }));
        assert_eq!(*active.borrow(), Some(marker.clone()));
        assert_eq!(*toy.calls.borrow(), vec![9.0 - WHEN_LOC_EPOCH_WINDOW_DAYS]);
        let toy = Toy {
            calls: RefCell::new(Vec::new()),
            global: marker.clone(),
            local: true,
        };
        assert_eq!(
            sol_eclipse_when_loc_source_via(&toy, 9.0, false, &|_| true, &|_| panic!(
                "local maximum errors are not global catches"
            )),
            Err(marker)
        );
    }
    #[test]
    fn raw_maximum_preserves_nan_comparison_and_terminal_provider_order() {
        let calls = RefCell::new(Vec::new());
        let mut first = true;
        let mut sample = |jd| {
            calls.borrow_mut().push(jd);
            if first {
                first = false;
                Ok(f64::NAN)
            } else {
                Err(Error::UnknownBody { body_id: 777 })
            }
        };
        assert_eq!(
            find_local_maximum_source_via(&mut sample, 10.0, 0.125),
            Err(Error::UnknownBody { body_id: 777 })
        );
        assert_eq!(calls.borrow().len(), 2);
        assert!(calls.borrow()[0] < calls.borrow()[1]);
    }
}
