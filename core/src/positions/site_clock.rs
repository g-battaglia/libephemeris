//! Configured topocentric clock pairing for one immutable public entry.

use crate::errors::{Error, Result};
use crate::time::{DeltaTConfig, DeltaTModel, deltat_with_config};

/// Pair a TT sample with the configured UT1 clock used by its site rotation.
///
/// All epochs and returned offsets are scalar Julian days. For a TT entry
/// (`original_ut == None`), with configured Delta-T in days `D`, evaluate
/// `d0 = D(T)` then return `T - D(T - d0)`. This is one fixed-point refinement,
/// not an inverse of the independently bundled frame or lunar timescale.
///
/// A UT1 entry carries its original scalar `U`. Form `T0 = U + D(U)`; the
/// central TT word returns `U` directly. A different TT sample uses
/// `T - (T0 - U)`, freezing the central subtraction rather than resampling D.
/// The downstream topocentric reducer also freezes its central TT-UT1 offset
/// across its internal SPEED stencil. Entry metadata is immutable, and the
/// callback is only consumed for an actual site rotation after site validation.
///
/// The enclosing position entry owns finite caller-JD validation; computed
/// IEEE values retain the downstream checked rotation/sampler policy.
/// Configured lookup errors propagate unchanged. No generic typed error is
/// treated as permission to substitute a stored Delta-T table: doing so needs
/// an independently established recoverable lookup classification.
pub(super) fn configured_site_ut1_partner(
    tt: f64,
    original_ut: Option<f64>,
    delta_t_days: &dyn Fn(f64) -> Result<f64>,
) -> Result<f64> {
    if let Some(ut) = original_ut {
        let central_tt = ut + delta_t_days(ut)?;
        if tt == central_tt {
            Ok(ut)
        } else {
            Ok(tt - (central_tt - ut))
        }
    } else {
        let initial = delta_t_days(tt)?;
        Ok(tt - delta_t_days(tt - initial)?)
    }
}

/// Pair a stored-state site sample with its configured UT1 clock and table.
///
/// Scalar TT and UT1 epochs and Delta-T are in Julian days. An original UT1
/// entry retains the immutable central word and frozen sample difference of
/// [`configured_site_ut1_partner`]; it never invokes the stored table.
/// A TT entry evaluates `D(T)` and then `D(T-D(T))`. A successful second
/// evaluation gives `U=T-D(T-D(T))`.
///
/// The narrowly recoverable calendar domain is an estimated UT1 of NaN under
/// Espenak-Meeus, disabled IERS and no user override, refused as
/// [`Error::NonFiniteArgument`] for that NaN argument. In this domain use
/// `U=T-table(T)` from the actual selected reader. Infinity, first-lookup
/// failures, unavailable models, other error variants and other model/backend
/// configurations retain their actual errors. A raw user override is returned
/// by both configured evaluations verbatim and cannot activate this branch.
///
/// `stored_delta_t_days` must belong to the same acquired immutable reader as
/// the position, and follows that provider's own table/tier selection policy.
/// Its result or failure is propagated without another substitution. No
/// mutable clock state or bundled frame/time model is introduced. The caller
/// owns public finite-JD and observer preflight; downstream SPEED sampling
/// freezes the returned central TT-minus-UT1 difference. Recovery for other
/// backend errors has no classification in this helper.
pub(crate) fn configured_site_ut1_partner_with_table(
    tt: f64,
    original_ut: Option<f64>,
    cfg: &DeltaTConfig<'_>,
    stored_delta_t_days: &dyn Fn(f64) -> Result<f64>,
) -> Result<f64> {
    let configured = |jd| deltat_with_config(jd, cfg);
    if original_ut.is_some() {
        return configured_site_ut1_partner(tt, original_ut, &configured);
    }
    tt_site_partner_with_fallback(tt, cfg, &configured, stored_delta_t_days)
}

/// Pair a differentiating sample with the clock of its public entry.
///
/// A TT entry uses the two-lookup inverse above. A UT1 entry instead inverts
/// its forward curve `T(U)=U+D(U)` around the original `U`; applying the TT
/// entry's distinct inverse could select another calendar piece. Fixed-point
/// refinement stops at the binary64 forward solution or a stable word. An
/// unresolvable clock after sixteen refinements is [`Error::MathDomain`].
/// Epochs and Delta-T are in Julian days; actual configured/table errors
/// propagate. The caller restricts samples to their central clock support.
pub(crate) fn configured_site_sample_ut1(
    tt: f64,
    original_ut: Option<f64>,
    cfg: &DeltaTConfig<'_>,
    stored_delta_t_days: &dyn Fn(f64) -> Result<f64>,
) -> Result<f64> {
    let Some(original) = original_ut else {
        return configured_site_ut1_partner_with_table(tt, None, cfg, stored_delta_t_days);
    };
    let central_tt = original + deltat_with_config(original, cfg)?;
    if tt == central_tt {
        return Ok(original);
    }
    let mut ut = original + (tt - central_tt);
    let mut previous = None;
    for _ in 0..16 {
        let delta = deltat_with_config(ut, cfg)?;
        if ut + delta == tt {
            return Ok(ut);
        }
        let next = tt - delta;
        if next == ut {
            return Ok(ut);
        }
        if previous == Some(next) {
            // Adjacent words can straddle the rounded inverse. Keep the one
            // with the smaller actual forward-clock residual.
            let next_error = (next + deltat_with_config(next, cfg)? - tt).abs();
            return Ok(if next_error < (ut + delta - tt).abs() {
                next
            } else {
                ut
            });
        }
        previous = Some(ut);
        ut = next;
    }
    Err(Error::MathDomain { value: tt })
}

/// Continuous TT offsets of the calendar pieces consumed by a site inverse.
///
/// The Espenak-Meeus selector uses the month's midpoint decimal year. Its
/// value consequently jumps at Gregorian month boundaries. Both `D(T)` and
/// `D(T-D(T))` must stay in their central months when differentiating a place.
/// For a UT1 entry, the support is instead the TT image of the original
/// forward-clock month. A user override or another clock model has no such
/// calendar restriction. Returned bounds are offsets in days, including their actual binary64 join
/// ownership. Bisection only locates a calendar boundary; it never reads a
/// position, speed or comparison residual.
pub(crate) fn configured_site_clock_support(
    tt: f64,
    original_ut: Option<f64>,
    cfg: &DeltaTConfig<'_>,
) -> Result<(f64, f64)> {
    if cfg.user_override_days.is_some()
        || cfg.model != DeltaTModel::EspenakMeeus
        || cfg.iers_enabled
        || !cfg.tidal_acceleration.is_finite()
    {
        return Ok((f64::NEG_INFINITY, f64::INFINITY));
    }
    use crate::time::{Calendar, julday, revjul};
    let calendar = revjul(original_ut.unwrap_or(tt), Calendar::Gregorian)?;
    let start = julday(calendar.year, calendar.month, 1, 0.0, Calendar::Gregorian)?;
    let (year, month) = if calendar.month == 12 {
        (
            calendar.year.checked_add(1).ok_or(Error::IntegerOverflow)?,
            1,
        )
    } else {
        (calendar.year, calendar.month + 1)
    };
    let end = julday(year, month, 1, 0.0, Calendar::Gregorian)?.next_down();
    if original_ut.is_some() {
        // Image of the one owned forward-clock month. The first endpoint is
        // included; the final representable UT word remains in that month.
        return Ok((
            start + deltat_with_config(start, cfg)? - tt,
            end + deltat_with_config(end, cfg)? - tt,
        ));
    }
    let refined = revjul(tt - deltat_with_config(tt, cfg)?, Calendar::Gregorian)?;
    let owned = |epoch| -> Result<bool> {
        let date = revjul(epoch - deltat_with_config(epoch, cfg)?, Calendar::Gregorian)?;
        Ok((date.year, date.month) == (refined.year, refined.month))
    };
    let boundary = |mut outside: f64, mut inside: f64| -> Result<f64> {
        while outside != inside {
            let middle = outside + (inside - outside) / 2.0;
            if middle == outside || middle == inside {
                break;
            }
            if owned(middle)? {
                inside = middle;
            } else {
                outside = middle;
            }
        }
        Ok(inside)
    };
    let lower = if owned(start)? {
        start
    } else {
        boundary(start, tt)?
    };
    let upper = if owned(end)? { end } else { boundary(end, tt)? };
    Ok((lower - tt, upper - tt))
}

/// Evaluate the two TT inverse lookups and their classified table alternative.
///
/// `configured_delta_t_days` realizes the immutable `cfg` priority chain;
/// the separate table callback serves its reader's Delta-T in days at original
/// TT. Only a second NaN-argument calendar refusal in the stated Espenak domain
/// is recoverable. The first error and every unclassified second error remain
/// owned and unchanged; a consumed table error replaces the caught failure.
fn tt_site_partner_with_fallback(
    tt: f64,
    cfg: &DeltaTConfig<'_>,
    configured_delta_t_days: &dyn Fn(f64) -> Result<f64>,
    stored_delta_t_days: &dyn Fn(f64) -> Result<f64>,
) -> Result<f64> {
    let initial = configured_delta_t_days(tt)?;
    let estimate = tt - initial;
    match configured_delta_t_days(estimate) {
        Ok(delta) => Ok(tt - delta),
        Err(error) => {
            let calendar_nan = cfg.model == DeltaTModel::EspenakMeeus
                && !cfg.iers_enabled
                && cfg.user_override_days.is_none()
                && estimate.is_nan()
                && matches!(&error, Error::NonFiniteArgument { value } if value.is_nan());
            if calendar_nan {
                Ok(tt - stored_delta_t_days(tt)?)
            } else {
                Err(error)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::errors::Error;
    use std::sync::Mutex;

    #[test]
    fn calendar_support_preserves_both_consumed_months_and_join_ownership() {
        use crate::time::{Calendar, julday, revjul};
        let cfg = DeltaTConfig {
            model: DeltaTModel::EspenakMeeus,
            iers_enabled: false,
            ..DeltaTConfig::default()
        };
        let tt = julday(2149, 12, 1, 0.0, Calendar::Gregorian).unwrap();
        let (lower, upper) = configured_site_clock_support(tt, None, &cfg).unwrap();
        assert_eq!(lower, 0.0);
        assert!(upper > 0.0 && upper < 1.0 / 64.0);
        let month_pair = |epoch: f64| {
            let first = revjul(epoch, Calendar::Gregorian).unwrap();
            let refined = revjul(
                epoch - deltat_with_config(epoch, &cfg).unwrap(),
                Calendar::Gregorian,
            )
            .unwrap();
            ((first.year, first.month), (refined.year, refined.month))
        };
        let central = month_pair(tt);
        assert_eq!(central, ((2149, 12), (2149, 11)));
        assert_eq!(month_pair(tt + lower), central);
        assert_eq!(month_pair(tt + upper), central);
        assert_ne!(month_pair(tt.next_down()), central);
        assert_ne!(month_pair((tt + upper).next_up()), central);
    }

    #[test]
    fn ut_samples_invert_their_own_forward_month_with_tidal_adjustment() {
        use crate::time::{Calendar, julday, revjul};
        let cfg = DeltaTConfig {
            model: DeltaTModel::EspenakMeeus,
            tidal_acceleration: -20.0,
            iers_enabled: false,
            ..DeltaTConfig::default()
        };
        let original = julday(2149, 11, 30, 23.99, Calendar::Gregorian).unwrap();
        let tt = original + deltat_with_config(original, &cfg).unwrap();
        let support = configured_site_clock_support(tt, Some(original), &cfg).unwrap();
        assert_eq!(revjul(tt, Calendar::Gregorian).unwrap().month, 12);
        for offset in [support.0, -1.0 / 1024.0, 0.0, support.1] {
            let sample = tt + offset;
            let ut = configured_site_sample_ut1(sample, Some(original), &cfg, &|_| {
                panic!("no stored table in the finite forward-clock inverse")
            })
            .unwrap();
            assert_eq!(revjul(ut, Calendar::Gregorian).unwrap().month, 11);
            assert_eq!(ut + deltat_with_config(ut, &cfg).unwrap(), sample);
            if offset == 0.0 {
                assert_eq!(ut.to_bits(), original.to_bits());
            }
        }
    }

    #[test]
    fn hard_clock_override_has_no_calendar_piece_boundary() {
        let cfg = DeltaTConfig {
            model: DeltaTModel::EspenakMeeus,
            user_override_days: Some(64.0 / 86400.0),
            ..DeltaTConfig::default()
        };
        assert_eq!(
            configured_site_clock_support(2506300.5, None, &cfg).unwrap(),
            (f64::NEG_INFINITY, f64::INFINITY)
        );
    }

    #[test]
    fn tt_refines_at_the_first_estimated_ut_before_final_subtraction() {
        let calls = Mutex::new(Vec::new());
        let d = |jd| {
            calls.lock().unwrap().push(jd);
            Ok(if jd == 100.0 { 0.25 } else { 0.5 })
        };
        assert_eq!(configured_site_ut1_partner(100.0, None, &d).unwrap(), 99.5);
        assert_eq!(*calls.lock().unwrap(), [100.0, 99.75]);
    }

    #[test]
    fn ut_retains_the_original_word_even_when_subtraction_cannot_recover_it() {
        let ut: f64 = 1.0;
        let delta: f64 = 1.0e16;
        let tt = ut + delta;
        assert_ne!((tt - delta).to_bits(), ut.to_bits());
        let d = |jd: f64| {
            assert_eq!(jd.to_bits(), ut.to_bits());
            Ok(delta)
        };
        assert_eq!(
            configured_site_ut1_partner(tt, Some(ut), &d)
                .unwrap()
                .to_bits(),
            ut.to_bits()
        );
    }

    #[test]
    fn ut_side_samples_freeze_the_actual_central_subtraction() {
        let ut: f64 = 2451545.125;
        let delta: f64 = 0.000_742_145;
        let tt = ut + delta;
        let frozen = tt - ut;
        let d = |jd: f64| {
            assert_eq!(jd.to_bits(), ut.to_bits());
            Ok(delta)
        };
        for sample in [tt - 0.0005, tt + 0.0005] {
            assert_eq!(
                configured_site_ut1_partner(sample, Some(ut), &d)
                    .unwrap()
                    .to_bits(),
                (sample - frozen).to_bits()
            );
        }
    }

    #[test]
    fn tt_first_lookup_error_is_terminal_without_second_lookup() {
        let calls = Mutex::new(0);
        let d = |_| {
            *calls.lock().unwrap() += 1;
            Err(Error::DeltaTModelUnavailable { model: "smh2016" })
        };
        assert!(matches!(
            configured_site_ut1_partner(2451545.0, None, &d),
            Err(Error::DeltaTModelUnavailable { model: "smh2016" })
        ));
        assert_eq!(*calls.lock().unwrap(), 1);
    }

    #[test]
    fn tt_second_lookup_error_preserves_its_actual_error() {
        let calls = Mutex::new(0);
        let d = |_| {
            let mut count = calls.lock().unwrap();
            *count += 1;
            if *count == 1 {
                Ok(0.25)
            } else {
                Err(Error::NonFiniteArgument {
                    value: f64::INFINITY,
                })
            }
        };
        assert!(matches!(
            configured_site_ut1_partner(100.0, None, &d),
            Err(Error::NonFiniteArgument { value }) if value == f64::INFINITY
        ));
        assert_eq!(*calls.lock().unwrap(), 2);
    }

    #[test]
    fn raw_user_override_ieee_values_reach_the_existing_sampler_boundary() {
        let d = |_| Ok(f64::NAN);
        assert!(
            configured_site_ut1_partner(2451545.0, None, &d)
                .unwrap()
                .is_nan()
        );
        assert!(
            configured_site_ut1_partner(2451545.0, Some(2451545.0), &d)
                .unwrap()
                .is_nan()
        );
    }
}

#[cfg(test)]
mod canonical_fallback_tests {
    use super::*;
    use std::cell::RefCell;

    /// Explicit finite-calendar model without optional provider acquisition.
    fn espenak() -> DeltaTConfig<'static> {
        DeltaTConfig {
            model: DeltaTModel::EspenakMeeus,
            iers_enabled: false,
            ..DeltaTConfig::default()
        }
    }

    #[test]
    fn computed_nan_inverse_uses_table_at_original_tt_and_source_partner_word() {
        let tt = f64::from_bits(0x4142492c10000000);
        let cfg = DeltaTConfig {
            tidal_acceleration: f64::NAN,
            ..espenak()
        };
        let first = deltat_with_config(tt, &cfg).unwrap();
        assert!(first.is_nan());
        assert!(matches!(deltat_with_config(tt - first, &cfg),
            Err(Error::NonFiniteArgument { value }) if value.is_nan()));
        let queries = RefCell::new(Vec::new());
        let table = |jd: f64| {
            queries.borrow_mut().push(jd.to_bits());
            Ok(f64::from_bits(0x3f13dff3975da36b))
        };
        let ut = configured_site_ut1_partner_with_table(tt, None, &cfg, &table).unwrap();
        assert_eq!(ut.to_bits(), 0x4142492c0ffd8402);
        assert_eq!(*queries.borrow(), [tt.to_bits()]);
        let frozen = tt - ut;
        for (sample, expected) in [
            (0x4142492c0fef9db2, 0x4142492c0fed21b4),
            (0x4142492c1010624e, 0x4142492c100de650),
        ] {
            assert_eq!((f64::from_bits(sample) - frozen).to_bits(), expected);
        }
    }

    #[test]
    fn ordinary_finite_inverse_does_not_query_table() {
        let calls = RefCell::new(Vec::new());
        let d = |jd: f64| {
            calls.borrow_mut().push(jd);
            Ok(if jd == 100.0 { 0.25 } else { 0.5 })
        };
        let result = tt_site_partner_with_fallback(100.0, &espenak(), &d, &|_| {
            panic!("unused table on successful canonical lookup")
        })
        .unwrap();
        assert_eq!(result, 99.5);
        assert_eq!(*calls.borrow(), [100.0, 99.75]);
    }

    #[test]
    fn first_refusal_is_terminal_without_second_or_table_read() {
        let calls = RefCell::new(0);
        let d = |_| {
            *calls.borrow_mut() += 1;
            Err(Error::DeltaTModelUnavailable {
                model: "retained first refusal",
            })
        };
        let result = tt_site_partner_with_fallback(100.0, &espenak(), &d, &|_| {
            panic!("unused table after unclassified first refusal")
        });
        assert!(matches!(
            result,
            Err(Error::DeltaTModelUnavailable {
                model: "retained first refusal"
            })
        ));
        assert_eq!(*calls.borrow(), 1);
    }

    #[test]
    fn nonfinite_error_for_a_finite_second_argument_is_not_calendar_nan() {
        let d = |jd: f64| {
            if jd == 100.0 {
                Ok(0.25)
            } else {
                Err(Error::NonFiniteArgument { value: f64::NAN })
            }
        };
        let result = tt_site_partner_with_fallback(100.0, &espenak(), &d, &|_| {
            panic!("a generic nonfinite error is not permission to query the table")
        });
        assert!(matches!(result, Err(Error::NonFiniteArgument { value }) if value.is_nan()));
    }

    #[test]
    fn computed_nan_does_not_reclassify_unavailable_or_integer_errors() {
        for error in [
            Error::DeltaTModelUnavailable {
                model: "retained second refusal",
            },
            Error::IntegerOverflow,
        ] {
            let d = |jd: f64| {
                if jd == 100.0 {
                    Ok(f64::NAN)
                } else {
                    Err(error.clone())
                }
            };
            let result = tt_site_partner_with_fallback(100.0, &espenak(), &d, &|_| {
                panic!("unclassified second error must be retained")
            });
            assert_eq!(result.unwrap_err(), error);
        }
    }

    #[test]
    fn infinity_calendar_domain_is_terminal() {
        let d = |jd: f64| {
            if jd == 100.0 {
                Ok(f64::INFINITY)
            } else {
                Err(Error::NonFiniteArgument { value: jd })
            }
        };
        let result = tt_site_partner_with_fallback(100.0, &espenak(), &d, &|_| {
            panic!("infinity is not a recoverable NaN calendar argument")
        });
        assert!(
            matches!(result, Err(Error::NonFiniteArgument { value }) if value == f64::NEG_INFINITY)
        );
    }

    #[test]
    fn raw_override_nan_and_infinity_do_not_activate_table() {
        for user in [f64::NAN, f64::INFINITY] {
            let cfg = DeltaTConfig {
                user_override_days: Some(user),
                ..espenak()
            };
            let result = configured_site_ut1_partner_with_table(2451545.125, None, &cfg, &|_| {
                panic!("override returned verbatim on both configured lookups")
            })
            .unwrap();
            if user.is_nan() {
                assert!(result.is_nan());
            } else {
                assert_eq!(result, f64::NEG_INFINITY);
            }
        }
    }

    #[test]
    fn original_ut_word_and_nan_path_never_use_tt_table() {
        let cfg = DeltaTConfig {
            user_override_days: Some(1e16),
            ..espenak()
        };
        let ut = 1.0;
        let tt = ut + 1e16;
        let value = configured_site_ut1_partner_with_table(tt, Some(ut), &cfg, &|_| {
            panic!("UT entry never consumes the TT table alternative")
        })
        .unwrap();
        assert_eq!(value.to_bits(), ut.to_bits());
        assert_ne!(value.to_bits(), (tt - 1e16).to_bits());
        let cfg = DeltaTConfig {
            tidal_acceleration: f64::NAN,
            ..espenak()
        };
        let value =
            configured_site_ut1_partner_with_table(2396760.125, Some(2396760.125), &cfg, &|_| {
                panic!("UT computed NaN does not enable the TT alternative")
            })
            .unwrap();
        assert!(value.is_nan());
    }

    #[test]
    fn other_models_and_enabled_iers_have_no_recoverable_classification() {
        for cfg in [
            DeltaTConfig {
                iers_enabled: true,
                ..espenak()
            },
            DeltaTConfig::default(),
        ] {
            let d = |jd: f64| {
                if jd == 100.0 {
                    Ok(f64::NAN)
                } else {
                    Err(Error::NonFiniteArgument { value: jd })
                }
            };
            let result = tt_site_partner_with_fallback(100.0, &cfg, &d, &|_| {
                panic!("unclassified model/backend configurations remain terminal")
            });
            assert!(matches!(result, Err(Error::NonFiniteArgument { value }) if value.is_nan()));
        }
    }

    #[test]
    fn actual_selected_table_failure_is_propagated() {
        let cfg = DeltaTConfig {
            tidal_acceleration: f64::NAN,
            ..espenak()
        };
        let error = Error::LebCorrupt {
            what: "selected table marker".into(),
            offset: 17,
        };
        let result =
            configured_site_ut1_partner_with_table(
                2396760.125,
                None,
                &cfg,
                &|_| Err(error.clone()),
            );
        assert_eq!(result.unwrap_err(), error);
    }
}
