//! Preparation, lazy delegation and ordered UT1 fixed-star batches.
use super::*;

/// Complete ordinary fast-batch lookup before its shared clock/observer.
/// Original keys and skipped slots are owned. Construction retains only the
/// actual terminal resolver failure; no provider or Session is acquired here.
#[derive(Debug)]
pub struct PreparedStarBatch {
    keys: Vec<String>,
    skip_errors: bool,
    any_resolved: bool,
}

/// Original geocentric pre-resolution stage; no added per-key queries.
/// This token is only used by the ordinary shared-Time branch, never to make
/// delegated scalar lookup eager. Empty/all-skipped input consumes no clock.
pub fn prepare_star_batch(
    keys: &[&str],
    skip_errors: bool,
) -> RetainedStarResult<PreparedStarBatch> {
    let mut any_resolved = false;
    for key in keys {
        match resolve_star_ref(key) {
            Ok(_) => any_resolved = true,
            Err(StarResolveError::Miss(_)) if skip_errors => {}
            Err(StarResolveError::Miss(failure)) => {
                return RetainedStarResult::failure(failure.to_error(key));
            }
            Err(StarResolveError::Refused(error)) => return RetainedStarResult::failure(error),
        }
    }
    RetainedStarResult::success(PreparedStarBatch {
        keys: keys.iter().map(|s| (*s).to_owned()).collect(),
        skip_errors,
        any_resolved,
    })
}

impl PreparedStarBatch {
    /// Empty/all-skipped result without context/configuration acquisition.
    pub fn skipped_slots(&self) -> Option<usize> {
        if self.any_resolved {
            None
        } else {
            Some(self.keys.len())
        }
    }
}

impl Session {
    /// Whether this batch uses lazy scalar delegation rather than shared Time.
    /// Reads only this Session's sidereal selection and input flags; no clock,
    /// provider or model/path setting is acquired.
    pub fn star_batch_delegates(&self, flags: i32) -> bool {
        let (selection, _) = crate::houses::sidereal_selection_for_context(&self.ctx);
        flags & FLG_TOPOCTR != 0
            || (flags & FLG_SIDEREAL != 0
                && (crate::sidereal::is_fixed_epoch_request(flags, selection.mode())
                    || (selection.bits() & (SIDBIT_ECL_T0 | SIDBIT_SSY_PLANE) != 0
                        && !crate::sidereal::SIDBIT_PROJECTION_SUPPRESS_MODES
                            .contains(&selection.mode()))))
    }

    /// V1 catalog places at one original UT1, preserving input order and slots.
    ///
    /// Each success is an owned [`FixedStar`] in the requested units/frame with
    /// UT normalized returned flags. With `skip_errors=false`, the first failure
    /// aborts. With skipping, catchable per-star library failures occupy `None`.
    /// Unicode decimal indices have no fixed digit-count limit; overflowing or
    /// out-of-range values are catalog misses, as are unresolved nondecimal
    /// name spellings. No partial vector is returned on terminal failure. Lookup follows [`Self::fixstar`], with no V2 wildcard.
    ///
    /// Topocentric, fixed-epoch sidereal and active projection requests resolve
    /// and evaluate each name lazily through [`Self::fixstar_ut`]. Their site,
    /// clock and provider failures belong to the current key and can be skipped;
    /// empty input consumes none of those dependencies, even for TOPOCTR.
    ///
    /// The ordinary geocentric route resolves all names before acquiring one
    /// shared observer at the original bundled UT1 Time, retaining its split TT
    /// and frame clock. If no names resolve, skipped slots return before time
    /// admission. SPEED shares fresh scalar-TT observers at TT +/- 0.5 days.
    /// Shared observer/clock failures remain terminal; only failures inside the
    /// per-star reduction loop can be skipped. A consumed raw shared-observer
    /// range error uses the original central TT envelope; per-star errors retain
    /// their raw identity. This retained-Time channel can differ in binary64
    /// words from repeated scalar UT calls and is preserved deliberately.
    /// Model Delta-T remains separate for sidereal defining epochs. LEB retains
    /// its stored per-key route; no provider or error-family substitute is added.
    pub fn batch_fixstars_ut(
        &self,
        keys: &[&str],
        jd_ut: f64,
        flags: i32,
        skip_errors: bool,
    ) -> Result<Vec<Option<FixedStar>>, Error> {
        self.batch_fixstars_ut_captured(keys, jd_ut, flags, skip_errors)
            .result
    }

    /// UT1 batch with its ordinary native terminal result.
    ///
    /// Calculation, slot order, units, clocks and skip policy are exactly those
    /// of [`Self::batch_fixstars_ut`]. Delegated keys are evaluated lazily;
    /// ordinary names resolve before the shared observer. Skipped errors and
    /// successful/empty batches carry no sidecar. Actual catalog/refusal fields
    /// are retained before conversion; acquired ranges use the serving filter.
    /// No partial result, backend exception object or fabricated graph escapes.
    pub fn batch_fixstars_ut_captured(
        &self,
        keys: &[&str],
        jd_ut: f64,
        flags: i32,
        skip_errors: bool,
    ) -> RetainedStarResult<Vec<Option<FixedStar>>> {
        let (selection, _) = crate::houses::sidereal_selection_for_context(&self.ctx);
        let delegate = flags & FLG_TOPOCTR != 0
            || (flags & FLG_SIDEREAL != 0
                && (crate::sidereal::is_fixed_epoch_request(flags, selection.mode())
                    || (selection.bits() & (SIDBIT_ECL_T0 | SIDBIT_SSY_PLANE) != 0
                        && !crate::sidereal::SIDBIT_PROJECTION_SUPPRESS_MODES
                            .contains(&selection.mode()))));
        if delegate {
            return scalar_batch_retained_via(keys, skip_errors, |key| {
                self.fixed_star_entry_with_role(key, jd_ut, flags, true, false)
            });
        }
        let prepared = prepare_star_batch(keys, skip_errors);
        match prepared.result {
            Ok(prepared) => self.batch_fixstars_prepared_ut_captured(prepared, jd_ut, flags),
            Err(error) => RetainedStarResult::failure(error),
        }
    }

    /// Ordinary shared-Time branch for an actual prepared batch. Original
    /// serving lookup, slot reduction, SPEED samples and terminal ranges remain
    /// unchanged. Delegated inputs must instead use lazy scalar calls; callers
    /// test `star_batch_delegates` over this same retained state first.
    pub fn batch_fixstars_prepared_ut_captured(
        &self,
        prepared: PreparedStarBatch,
        jd_ut: f64,
        flags: i32,
    ) -> RetainedStarResult<Vec<Option<FixedStar>>> {
        let keys: Vec<&str> = prepared.keys.iter().map(String::as_str).collect();
        let keys = keys.as_slice();
        let skip_errors = prepared.skip_errors;
        if !prepared.any_resolved {
            return RetainedStarResult::success(vec![None; keys.len()]);
        }
        let selected = match crate::positions::star_serving::auto_star_context(
            &self.ctx, jd_ut, true, flags,
        ) {
            Ok(selected) => selected,
            Err(error) => return RetainedStarResult::failure(error),
        };
        let context = selected.as_ref().unwrap_or(&self.ctx);
        let acquired = crate::positions::star_serving::with_context_star_serving_captured_with_main(
            context,
            &self.star_main,
            None,
            self.delta_t_model(),
            self.ctx.effective_iers_delta_t(),
            flags,
            self.ctx
                .observer()
                .map(|site| (site.lon_deg(), site.lat_deg(), site.height_m())),
            |serving, topo| {
                serving
                    .batch_fixstars_ut(keys, jd_ut, flags, topo, skip_errors)
                    .map(|rows| {
                        rows.into_iter()
                            .enumerate()
                            .map(|(index, row)| {
                                row.map(|(values, canonical_name, returned_flags)| {
                                    if let Ok(star) = resolve_star_ref(keys[index]) {
                                        crate::trace::record(star.star_id, serving.source_tag());
                                    }
                                    FixedStar {
                                        values,
                                        canonical_name,
                                        returned_flags,
                                    }
                                })
                            })
                            .collect()
                    })
            },
        );
        RetainedStarResult::from_serving(acquired)
    }
}

/// Evaluate actual delegated scalar calls in input order and retain only a
/// terminal error. Empty inputs acquire nothing. With skipping, library errors
/// become None and their sidecars are dropped; converter refusals always abort.
/// Success returns the unchanged owned slots with no stale failure metadata.
fn scalar_batch_retained_via<T>(
    keys: &[&str],
    skip_errors: bool,
    mut scalar: impl FnMut(&str) -> Result<T, RetainedScalarStarFailure>,
) -> RetainedStarResult<Vec<Option<T>>> {
    let mut out = Vec::with_capacity(keys.len());
    for key in keys {
        match scalar(key) {
            Ok(value) => out.push(Some(value)),
            Err(failure)
                if skip_errors && matches!(&failure.role, ScalarStarFailure::Library(_)) =>
            {
                out.push(None)
            }
            Err(failure) => return failure.into_public(),
        }
    }
    RetainedStarResult::success(out)
}

/// Evaluate a delegated batch lazily in input order. The callback owns catalog,
/// configuration and scalar clock acquisition. Empty input never invokes it;
/// skipped library failures occupy `None`, while integer-conversion refusals
/// abort even with skipping. Successful owned results escape only on completion.
#[cfg(test)]
fn scalar_batch_via<T>(
    keys: &[&str],
    skip_errors: bool,
    mut scalar: impl FnMut(&str) -> Result<T, ScalarStarFailure>,
) -> Result<Vec<Option<T>>, Error> {
    scalar_batch_retained_via(keys, skip_errors, |key| {
        scalar(key).map_err(|role| RetainedScalarStarFailure { role })
    })
    .result
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn empty_topocentric_batch_consumes_no_site_clock_or_provider() {
        let session =
            Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
        for skip in [false, true] {
            assert_eq!(
                session.batch_fixstars_ut(&[], f64::NAN, FLG_TOPOCTR, skip),
                Ok(vec![])
            );
        }
    }
    #[test]
    fn delegated_missing_names_and_configuration_are_per_key() {
        let session =
            Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
        assert_eq!(
            session.batch_fixstars_ut(
                &["no-such-star-xyz", "Regulus"],
                2451545.125,
                FLG_TOPOCTR,
                true
            ),
            Ok(vec![None, None])
        );
        assert!(matches!(
            session.batch_fixstars_ut(
                &["no-such-star-xyz", "Regulus"],
                2451545.125,
                FLG_TOPOCTR,
                false
            ),
            Err(Error::StarNotFound { .. })
        ));
    }
    #[test]
    fn shared_skipped_misses_acquire_no_provider_or_clock() {
        let session =
            Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
        assert_eq!(
            session.batch_fixstars_ut(&["no-such-star-xyz"], f64::NAN, 0, true),
            Ok(vec![None])
        );
        assert_eq!(
            session.batch_fixstars_ut(&["²"], 2451545.125, 0, true),
            Ok(vec![None])
        );
    }
    #[test]
    fn delegated_callback_is_lazy_and_preserves_failure_order() {
        let mut calls = vec![];
        let result: Result<Vec<Option<&str>>, Error> =
            scalar_batch_via(&["missing", "configured", "later"], false, |key| {
                calls.push(key.to_owned());
                Err(ScalarStarFailure::library(Error::ConfigurationMissing {
                    missing_config: "observer_location",
                    suggestion: "synthetic refusing callback",
                }))
            });
        assert!(result.is_err());
        assert_eq!(calls, ["missing"]);
    }
    #[test]
    fn delegated_callback_keeps_successes_and_skips_library_failures() {
        let mut calls = vec![];
        let result = scalar_batch_via(&["missing", "ok", "configuration"], true, |key| {
            calls.push(key.to_owned());
            if key == "ok" {
                Ok(7)
            } else {
                Err(ScalarStarFailure::library(Error::NonFiniteArgument {
                    value: f64::NAN,
                }))
            }
        })
        .unwrap();
        assert_eq!(result, [None, Some(7), None]);
        assert_eq!(calls, ["missing", "ok", "configuration"]);
    }
    #[test]
    fn delegated_empty_callback_is_never_called() {
        let result = scalar_batch_via::<i32>(&[], false, |_| {
            panic!("empty delegated batch consumed callback")
        });
        assert_eq!(result, Ok(vec![]));
    }

    #[test]
    fn delegated_nondecimal_numeric_names_follow_missing_name_policy() {
        let session =
            Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
        assert!(matches!(
            session.batch_fixstars_ut(&["²", "later"], f64::NAN, FLG_TOPOCTR, false),
            Err(Error::StarNotFound { star_id, search_type: Some("name") }) if star_id == "²"
        ));
        assert_eq!(
            session.batch_fixstars_ut(&["²", "later"], f64::NAN, FLG_TOPOCTR, true),
            Ok(vec![None, None])
        );
    }

    #[test]
    fn delegated_late_refusal_aborts_after_current_key_without_prepass() {
        let mut calls = vec![];
        let result = scalar_batch_via(&["missing", "ok", "refused", "later"], true, |key| {
            calls.push(key.to_owned());
            match key {
                "missing" => Err(ScalarStarFailure::library(Error::StarNotFound {
                    star_id: "synthetic missing".to_owned(),
                    search_type: Some("name"),
                })),
                "refused" => Err(ScalarStarFailure::refused(Error::InvalidStarSearch {
                    input: key.to_owned(),
                })),
                _ => Ok(7),
            }
        });
        assert!(matches!(result, Err(Error::InvalidStarSearch { .. })));
        assert_eq!(calls, ["missing", "ok", "refused"]);
    }
}

#[cfg(test)]
mod prepared_star_boundary_tests {
    use super::*;

    #[test]
    fn all_skipped_preparation_never_consumes_session_or_epoch() {
        let prepared = prepare_star_batch(&["no-such-star-prepared-control"], true)
            .result
            .unwrap();
        assert_eq!(prepared.skipped_slots(), Some(1));
        let session = Session::from_context(
            crate::context::EphemerisContext::new().with_mode(ProviderMode::NativeJpl),
        );
        assert_eq!(
            session
                .batch_fixstars_prepared_ut_captured(prepared, f64::NAN, 0)
                .result
                .unwrap(),
            vec![None]
        );
    }

    #[test]
    fn nondecimal_numeric_name_is_a_skippable_catalog_miss() {
        let result = prepare_star_batch(&["no-such-star-prepared-control", "²"], true);
        assert_eq!(result.result.unwrap().skipped_slots(), Some(2));
    }
}
