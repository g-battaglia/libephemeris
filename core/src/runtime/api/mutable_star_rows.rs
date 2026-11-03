//! Fixed-star calculations over resolved, owned catalogue rows.
//!
//! Catalogue identity is retained separately from the per-sample numerical
//! row. Native provider, clock, correction and output-frame calculations
//! consume checked [`NativeStarInputs`]. Missing rows propagate acquisition
//! errors rather than selecting a static substitute.

use super::{FixedStar, RetainedStarResult, Session};
pub use crate::stars::{AcquiredStarRow, NativeStarInputs, StarRowAcquisition, StarRowBackend};
use std::sync::Arc;

/// Owned routing identity supplied by an actual caller catalog acquisition.
/// The canonical string is metadata, not a second resolver input. An attached
/// caller can keep the untouched raw name/entry separately for final formatting.
/// Construction does not read a row, epoch, provider, observer or configuration.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AcquiredFixedStar {
    hit: crate::stars::StarMatch,
    v2: bool,
}

impl AcquiredFixedStar {
    /// Retain the resolved numerical id and name, with the selected V1/V2 family.
    /// IDs outside the immutable catalog remain valid routing identities when
    /// the supplied live acquisition owns their actual row. No row is invented.
    pub fn new(star_id: i32, canonical_name: String, v2: bool) -> Self {
        Self {
            hit: crate::stars::StarMatch {
                star_id,
                canonical: canonical_name,
            },
            v2,
        }
    }

    /// Numerical routing id; the row adapter retains the original raw key.
    pub fn star_id(&self) -> i32 {
        self.hit.star_id
    }
}

impl Session {
    /// Scalar TT/UT1 place for a genuinely resolved live catalog identity.
    ///
    /// Input/output units, sidereal/coordinate flags and clocks follow the
    /// corresponding prepared fixed-star entry. Stored samples select their
    /// row before observer reads; native samples select after original Time and
    /// observer acquisition. SPEED first consumes membership-only velocity
    /// admission and then obtains separate central/minus/plus leases. Actual
    /// callback/provider failures are returned unchanged with owning range data.
    /// This method does not classify or render an attached exception.
    pub fn fixed_star_acquired_captured(
        &self,
        star: AcquiredFixedStar,
        jd: f64,
        flags: i32,
        is_ut: bool,
        rows: &dyn StarRowAcquisition,
    ) -> RetainedStarResult<FixedStar> {
        self.fixed_star_rows_captured(star.hit, jd, flags, is_ut, star.v2, rows)
    }

    /// Scalar fixed-star place with an optional live sidereal anchor catalog.
    ///
    /// Units, clocks, flags and acquired fixed-star leases follow
    /// [`Self::fixed_star_acquired_captured`]. The separate anchor supplier is
    /// borrowed only during this call. The fixed-star flag layer invokes it at
    /// the released `fixed_stars._positional_transform` nested ayanamsha stage,
    /// after geometry and selected frame operations. Center and SPEED endpoint
    /// reductions make fresh helper/node requests; no anchor row is frozen.
    /// Non-sidereal and data-independent modes never request its rows.
    /// `None` retains standalone static behavior; supplied failures propagate.
    pub fn fixed_star_acquired_captured_with_anchor_rows(
        &self,
        star: AcquiredFixedStar,
        jd: f64,
        flags: i32,
        is_ut: bool,
        rows: &dyn StarRowAcquisition,
        anchors: Option<&dyn crate::stars::AnchorRowAcquisition>,
    ) -> RetainedStarResult<FixedStar> {
        self.fixed_star_rows_captured_with_anchors(
            star.hit, jd, flags, is_ut, star.v2, rows, anchors,
        )
    }

    /// Scalar acquired place with its actual per-observation Time/Main owner.
    /// Each native scalar TT sample acquires before observer evaluation; stored
    /// samples consume no native owner. Returned resources own immutable bytes,
    /// and errors keep the existing row/catch/presentation boundaries.
    #[allow(clippy::too_many_arguments)]
    pub fn fixed_star_acquired_captured_with_anchor_rows_and_acquisition(
        &self,
        star: AcquiredFixedStar,
        jd: f64,
        flags: i32,
        is_ut: bool,
        rows: &dyn StarRowAcquisition,
        anchors: Option<&dyn crate::stars::AnchorRowAcquisition>,
        acquisition: Option<&crate::positions::star_serving::StarObservationAcquisition<'_>>,
    ) -> RetainedStarResult<FixedStar> {
        self.fixed_star_rows_captured_with_anchors_and_acquisition(
            star.hit,
            jd,
            flags,
            is_ut,
            star.v2,
            rows,
            anchors,
            acquisition,
        )
    }

    /// Shared UT1 batch for already acquired V1 catalog identities.
    ///
    /// None slots are genuine prior skipped resolver misses. With no acquired
    /// star this returns before provider/Time selection. Native shared observers
    /// are acquired before any row; each star then gets independent central and
    /// speed-side leases without scalar velocity admission. Its actual per-star
    /// catch decision is fallible. Caught stored failures invoke their actual
    /// logging callback and retry the same retained LEB vector provider, keeping
    /// original Time and prior stored slots. Scalar delegated branches preserve
    /// their entry catches. No static name resolution or copied catalog occurs.
    /// Canonical result strings are retained metadata; raw name formatting,
    /// wider caller operators and backend/factory protocols belong to adapters.
    pub fn batch_fixstars_acquired_ut_captured(
        &self,
        resolved: &[Option<AcquiredFixedStar>],
        jd_ut: f64,
        flags: i32,
        skip_errors: bool,
        rows: &[Option<Arc<dyn StarRowAcquisition>>],
    ) -> RetainedStarResult<Vec<Option<FixedStar>>> {
        self.batch_fixstars_acquired_ut_captured_with_anchor_rows(
            resolved,
            jd_ut,
            flags,
            skip_errors,
            rows,
            None,
        )
    }

    /// Shared or delegated UT1 batch with independent live anchor acquisition.
    ///
    /// Entry units, occurrence alignment, shared observer/Time, retained partial
    /// stored results and catches follow
    /// [`Self::batch_fixstars_acquired_ut_captured`]. `anchors` remains borrowed
    /// through synchronous work and is requested only at each successful
    /// position's sidereal flag reduction. Delegated calls receive the same
    /// acquisition owner but still request fresh rows at their actual stages;
    /// the shared path preserves one observer and independent fixed-star
    /// occurrence owners. Empty/all-skipped inputs consume no anchor catalog.
    /// `None` retains the existing static catalog and error behavior.
    pub fn batch_fixstars_acquired_ut_captured_with_anchor_rows(
        &self,
        resolved: &[Option<AcquiredFixedStar>],
        jd_ut: f64,
        flags: i32,
        skip_errors: bool,
        rows: &[Option<Arc<dyn StarRowAcquisition>>],
        anchors: Option<&dyn crate::stars::AnchorRowAcquisition>,
    ) -> RetainedStarResult<Vec<Option<FixedStar>>> {
        self.batch_fixstars_acquired_ut_captured_with_anchor_rows_and_acquisition(
            resolved,
            jd_ut,
            flags,
            skip_errors,
            rows,
            anchors,
            None,
        )
    }

    /// Acquired batch with one real shared UT1 Time/Main acquisition, or fresh
    /// scalar-TT acquisitions in a genuinely delegated request. Empty/skipped
    /// batches consume none; source clocks and numerical reductions are intact.
    #[allow(clippy::too_many_arguments)]
    pub fn batch_fixstars_acquired_ut_captured_with_anchor_rows_and_acquisition(
        &self,
        resolved: &[Option<AcquiredFixedStar>],
        jd_ut: f64,
        flags: i32,
        skip_errors: bool,
        rows: &[Option<Arc<dyn StarRowAcquisition>>],
        anchors: Option<&dyn crate::stars::AnchorRowAcquisition>,
        acquisition: Option<&crate::positions::star_serving::StarObservationAcquisition<'_>>,
    ) -> RetainedStarResult<Vec<Option<FixedStar>>> {
        if rows.len() != resolved.len()
            || resolved
                .iter()
                .zip(rows)
                .any(|(star, row)| star.is_some() != row.is_some())
        {
            return RetainedStarResult::failure(super::Error::ConfigurationMissing {
                missing_config: "star_row_occurrence",
                suggestion: "each acquired batch identity requires its corresponding owned row acquisition",
            });
        }
        if resolved.iter().all(Option::is_none) {
            return RetainedStarResult::success(vec![None; resolved.len()]);
        }
        if self.star_batch_delegates(flags) {
            let mut output = Vec::with_capacity(resolved.len());
            for (star, rows) in resolved.iter().zip(rows) {
                let Some(star) = star else {
                    output.push(None);
                    continue;
                };
                let rows = rows.as_ref().expect("validated occurrence owner");
                let retained = self.fixed_star_acquired_captured_with_anchor_rows_and_acquisition(
                    star.clone(),
                    jd_ut,
                    flags,
                    true,
                    rows.as_ref(),
                    anchors,
                    acquisition,
                );
                match retained.result {
                    Ok(value) => output.push(Some(value)),
                    Err(error) => match rows.is_skippable(&error) {
                        Ok(true) if skip_errors => output.push(None),
                        Ok(_) => {
                            return RetainedStarResult::failure(error);
                        }
                        Err(classification) => {
                            return RetainedStarResult::failure(classification);
                        }
                    },
                }
            }
            return RetainedStarResult::success(output);
        }
        let matches = resolved
            .iter()
            .map(|star| star.as_ref().map(|star| star.hit.clone()))
            .collect::<Vec<_>>();
        let acquired =
            crate::positions::star_serving::with_context_star_serving_captured_with_rows_and_acquisition(
                &self.ctx,
                &self.star_main,
                None,
                self.delta_t_model(),
                self.ctx.effective_iers_delta_t(),
                flags,
                self.ctx
                    .observer()
                    .map(|site| (site.lon_deg(), site.lat_deg(), site.height_m())),
                None,
                acquisition,
                |serving, topo| {
                    let run = |provider: Option<&dyn crate::sidereal::DynamicAnchorProvider>| {
                        let live =
                            provider.map(|provider| serving.borrowing_dynamic_anchors(provider));
                        let serving = live.as_ref().unwrap_or(serving);
                        serving
                            .batch_resolved_ut(&matches, rows, jd_ut, flags, topo, skip_errors)
                            .map(|entries| {
                                entries
                                    .into_iter()
                                    .map(|entry| {
                                        entry.map(|(values, canonical_name, returned_flags)| {
                                            FixedStar {
                                                values,
                                                canonical_name,
                                                returned_flags,
                                            }
                                        })
                                    })
                                    .collect()
                            })
                    };
                    match anchors {
                        None => run(None),
                        Some(anchors) => {
                            let captured = self
                                .star_anchor_provider_captured(anchors, |provider| {
                                    run(Some(provider))
                                });
                            captured.result
                        }
                    }
                },
            );

        RetainedStarResult::from_serving(acquired)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::errors::Error;
    use crate::stars::{
        AcquiredAnchorRow, AnchorRowAcquisition, AnchorSample, AnchorSampleOptions,
        NativeStarInputs, SealedStarServing,
    };
    use std::sync::Mutex;

    struct UnusedProvider;
    impl crate::positions::StateProvider for UnusedProvider {
        fn provider_name(&self) -> &'static str {
            "unused cached-anchor states"
        }
        fn has_body(&self, _: i32) -> bool {
            panic!("completed cached anchor must not read body membership")
        }
        fn body_coverage(&self, _: i32) -> Option<(f64, f64)> {
            panic!("completed cached anchor must not read coverage")
        }
        fn serving_coord_type(&self, _: i32) -> Option<u32> {
            panic!("completed cached anchor must not read channel metadata")
        }
        fn eval_raw(&self, _: i32, _: f64) -> Result<(Vec<f64>, Vec<f64>), Error> {
            panic!("completed cached anchor must not read state")
        }
        fn eval_raw_two_part(
            &self,
            _: i32,
            _: crate::types::TwoPartEpoch,
        ) -> Result<(Vec<f64>, Vec<f64>), Error> {
            panic!("completed cached anchor must not read split state")
        }
        fn serving_entry_snapshot(&self, _: i32) -> Result<crate::positions::EntrySnapshot, Error> {
            panic!("completed cached anchor must not acquire an entry")
        }
    }

    struct CachedAnchorRows(Arc<Mutex<Vec<f64>>>);
    impl AnchorRowAcquisition for CachedAnchorRows {
        fn acquire_row(&self, key: &str) -> Result<Arc<dyn AcquiredAnchorRow>, Error> {
            assert_eq!(key, "SPICA");
            Ok(Arc::new(CachedAnchorRow(Arc::clone(&self.0))))
        }
    }

    struct CachedAnchorRow(Arc<Mutex<Vec<f64>>>);
    impl AcquiredAnchorRow for CachedAnchorRow {
        fn sample(&self, tt: f64, options: AnchorSampleOptions) -> Result<AnchorSample, Error> {
            assert!(options.nonut);
            self.0.lock().unwrap().push(tt);
            Ok(AnchorSample::Cached(195.0 + (tt - 2451545.0)))
        }
        fn galactic_inputs(&self) -> Result<NativeStarInputs, Error> {
            panic!("Spica shaping must not request Galactic Center operands")
        }
    }

    #[test]
    fn borrowed_live_anchor_shaping_reacquires_center_and_rate_endpoints() {
        use crate::constants::{FLG_J2000, FLG_SIDEREAL, FLG_SPEED, SIDM_TRUE_CITRA};
        let provider = UnusedProvider;
        let nutation = |_| panic!("native J2000 shape must not acquire nutation");
        let delta = |_| Ok(0.0);
        let mut selection = crate::sidereal::SiderealSelection::unset();
        selection.set_sidereal_mode(SIDM_TRUE_CITRA, 0.0, 0.0);
        let serving = SealedStarServing::new(&provider, &nutation, &delta, selection);
        let samples = Arc::new(Mutex::new(Vec::new()));
        let rows = CachedAnchorRows(Arc::clone(&samples));
        let live = serving.borrowing_live_anchor_rows(&rows);
        assert!(samples.lock().unwrap().is_empty());
        let shaped = live
            .shape_native_place(
                [100.0, 0.0, 1.0, 2.0, 0.0, 0.0],
                2451545.0,
                FLG_J2000 | FLG_SIDEREAL | FLG_SPEED,
                true,
            )
            .unwrap();
        // Exact input algebra: Citra offset is cached Spica longitude -180;
        // its one degree/day change subtracts from the raw two degree/day rate.
        assert_eq!(shaped, [85.0, 0.0, 1.0, 1.0, 0.0, 0.0]);
        assert_eq!(*samples.lock().unwrap(), [2451545.0, 2451545.5, 2451544.5]);
    }

    #[test]
    fn live_anchor_refusal_is_terminal_after_earlier_shape_admission() {
        use crate::constants::{FLG_J2000, FLG_SIDEREAL, SIDM_TRUE_CITRA};
        struct Refused;
        impl AnchorRowAcquisition for Refused {
            fn acquire_row(&self, _: &str) -> Result<Arc<dyn AcquiredAnchorRow>, Error> {
                Err(Error::ConfigurationMissing {
                    missing_config: "actual anchor callback refusal",
                    suggestion: "retain this owned refusal",
                })
            }
        }
        let provider = UnusedProvider;
        let nutation = |_| panic!("native J2000 shape must not acquire nutation");
        let delta = |_| Ok(0.0);
        let mut selection = crate::sidereal::SiderealSelection::unset();
        selection.set_sidereal_mode(SIDM_TRUE_CITRA, 0.0, 0.0);
        let serving = SealedStarServing::new(&provider, &nutation, &delta, selection);
        let rows = Refused;
        let live = serving.borrowing_live_anchor_rows(&rows);
        let flags = FLG_J2000 | FLG_SIDEREAL;
        assert!(matches!(
            live.shape_native_place([f64::NAN, 0.0, 1.0, 0.0, 0.0, 0.0], 2451545.0, flags, true),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            live.shape_native_place([100.0, 0.0, 1.0, 0.0, 0.0, 0.0], 2451545.0, flags, true),
            Err(Error::ConfigurationMissing {
                missing_config: "actual anchor callback refusal",
                ..
            })
        ));
    }

    #[test]
    fn independent_auto_and_stored_cached_anchors_need_no_main_or_reader() {
        use crate::constants::{FLG_J2000, FLG_SIDEREAL, SIDM_TRUE_CITRA};
        for mode in [
            crate::context::ProviderMode::Auto,
            crate::context::ProviderMode::SealedLeb,
        ] {
            let session = Session::from_context(
                crate::context::EphemerisContext::new()
                    .with_mode(mode)
                    .with_leb_file(None),
            );
            let provider = UnusedProvider;
            let nutation = |_| panic!("native J2000 shape must not acquire nutation");
            let delta = |_| Ok(0.0);
            let mut selection = crate::sidereal::SiderealSelection::unset();
            selection.set_sidereal_mode(SIDM_TRUE_CITRA, 0.0, 0.0);
            let body = SealedStarServing::new(&provider, &nutation, &delta, selection);
            let samples = Arc::new(Mutex::new(Vec::new()));
            let rows = CachedAnchorRows(Arc::clone(&samples));
            let result = session.star_anchor_provider_captured(&rows, |anchors| {
                body.borrowing_dynamic_anchors(anchors).shape_native_place(
                    [100.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                    2451545.0,
                    FLG_J2000 | FLG_SIDEREAL,
                    true,
                )
            });
            assert_eq!(result.result.unwrap(), [85.0, 0.0, 1.0, 0.0, 0.0, 0.0]);
            assert_eq!(*samples.lock().unwrap(), [2451545.0]);
        }
    }

    #[test]
    fn resolved_native_callback_receives_the_current_independent_anchor_owner() {
        use crate::constants::{FLG_J2000, FLG_SIDEREAL, SIDM_TRUE_CITRA};
        let session = Session::from_context(
            crate::context::EphemerisContext::new()
                .with_mode(crate::context::ProviderMode::Auto)
                .with_leb_file(None),
        );
        let provider = UnusedProvider;
        let nutation = |_| panic!("native J2000 shape must not acquire nutation");
        let delta = |_| Ok(0.0);
        let mut selection = crate::sidereal::SiderealSelection::unset();
        selection.set_sidereal_mode(SIDM_TRUE_CITRA, 0.0, 0.0);
        let body = SealedStarServing::new(&provider, &nutation, &delta, selection);
        let callback =
            |_: &[Option<crate::stars::StarMatch>],
             _: &[Option<Arc<dyn StarRowAcquisition>>],
             tt,
             flags,
             _,
             _,
             anchors: Option<&dyn crate::sidereal::DynamicAnchorProvider>| {
                let anchors = anchors.expect("resolved batch lost its live computation owner");
                let shaped = body.borrowing_dynamic_anchors(anchors).shape_native_place(
                    [100.0, 0.0, 1.0, 0.0, 0.0, 0.0],
                    tt,
                    flags,
                    true,
                )?;
                Ok(vec![Some((shaped, "actual occurrence".into(), flags))])
            };
        let batch = SealedStarServing::new(&provider, &nutation, &delta, selection)
            .with_resolved_native_anchor_batch(&callback);
        let samples = Arc::new(Mutex::new(Vec::new()));
        let rows = CachedAnchorRows(Arc::clone(&samples));
        let occurrences = [Some(crate::stars::StarMatch {
            star_id: 90000001,
            canonical: "actual occurrence".into(),
        })];
        let owners = [Some(Arc::new(Unused) as Arc<dyn StarRowAcquisition>)];
        let captured = session.star_anchor_provider_captured(&rows, |anchors| {
            batch.borrowing_dynamic_anchors(anchors).batch_resolved_ut(
                &occurrences,
                &owners,
                2451545.0,
                FLG_J2000 | FLG_SIDEREAL,
                None,
                false,
            )
        });
        assert_eq!(captured.result.unwrap()[0].as_ref().unwrap().0[0], 85.0);
        assert_eq!(*samples.lock().unwrap(), [2451545.0]);
    }

    struct Unused;
    impl AnchorRowAcquisition for Unused {
        fn acquire_row(&self, _: &str) -> Result<Arc<dyn AcquiredAnchorRow>, Error> {
            panic!("empty batch must not consume its live anchor supplier")
        }
    }
    impl StarRowAcquisition for Unused {
        fn velocity_admission(&self, _: i32) -> Result<(), super::super::Error> {
            panic!("unused")
        }
        fn acquire_row(
            &self,
            _: i32,
            _: StarRowBackend,
        ) -> Result<Arc<dyn AcquiredStarRow>, super::super::Error> {
            panic!("unused")
        }
        fn is_skippable(&self, _: &super::super::Error) -> Result<bool, super::super::Error> {
            panic!("unused")
        }
    }

    #[test]
    fn actual_observation_acquisition_precedes_native_observers_and_skips_empty_batches() {
        use crate::positions::star_serving::StarObservationClock;
        let session = Session::from_context(
            crate::context::EphemerisContext::new()
                .with_mode(crate::context::ProviderMode::NativeJpl),
        );
        let calls = Mutex::new(Vec::new());
        let marker = Error::NetworkSealed {
            purpose: "actual observer-stage factory refusal".into(),
        };
        let acquire = |clock| {
            calls.lock().unwrap().push(clock);
            Err(marker.clone())
        };
        let empty = session.batch_fixstars_acquired_ut_captured_with_anchor_rows_and_acquisition(
            &[None, None],
            f64::NAN,
            0,
            true,
            &[None, None],
            None,
            Some(&acquire),
        );
        assert_eq!(empty.result.unwrap(), vec![None, None]);
        assert!(calls.lock().unwrap().is_empty());
        let star = AcquiredFixedStar::new(90000001, String::new(), false);
        let scalar = session.fixed_star_acquired_captured_with_anchor_rows_and_acquisition(
            star.clone(),
            2451545.0,
            0,
            false,
            &Unused,
            None,
            Some(&acquire),
        );
        assert_eq!(scalar.result.unwrap_err().primary_error(), &marker);
        assert_eq!(
            *calls.lock().unwrap(),
            [StarObservationClock::Tt(2451545.0)]
        );
        calls.lock().unwrap().clear();
        let owners = [Some(Arc::new(Unused) as Arc<dyn StarRowAcquisition>)];
        let batch = session.batch_fixstars_acquired_ut_captured_with_anchor_rows_and_acquisition(
            &[Some(star)],
            2451545.0,
            0,
            true,
            &owners,
            None,
            Some(&acquire),
        );
        assert_eq!(batch.result.unwrap_err().primary_error(), &marker);
        assert_eq!(
            *calls.lock().unwrap(),
            [StarObservationClock::Ut1(2451545.0)]
        );
    }

    #[test]
    fn acquired_batch_all_misses_skips_nonfinite_clock_and_empty_main() {
        let session = Session::from_context(
            crate::context::EphemerisContext::new()
                .with_mode(crate::context::ProviderMode::NativeJpl),
        );
        let output = session.batch_fixstars_acquired_ut_captured(
            &[None, None],
            f64::NAN,
            crate::constants::FLG_TOPOCTR,
            true,
            &[None, None],
        );
        assert_eq!(output.result.unwrap(), vec![None, None]);
    }

    #[test]
    fn live_anchor_batch_all_misses_skips_catalog_clock_and_empty_main() {
        let session = Session::from_context(
            crate::context::EphemerisContext::new()
                .with_mode(crate::context::ProviderMode::NativeJpl),
        );
        let rows = Unused;
        let output = session.batch_fixstars_acquired_ut_captured_with_anchor_rows(
            &[None, None],
            f64::NAN,
            crate::constants::FLG_SIDEREAL | crate::constants::FLG_TOPOCTR,
            true,
            &[None, None],
            Some(&rows),
        );
        assert_eq!(output.result.unwrap(), vec![None, None]);
    }

    #[test]
    fn occurrence_alignment_refuses_before_any_callback_or_clock() {
        let session = Session::from_context(
            crate::context::EphemerisContext::new()
                .with_mode(crate::context::ProviderMode::NativeJpl),
        );
        let actual = AcquiredFixedStar::new(90000001, "actual acquired".into(), false);
        let output = session.batch_fixstars_acquired_ut_captured(
            &[Some(actual.clone())],
            f64::NAN,
            0,
            false,
            &[],
        );
        assert!(matches!(
            output.result,
            Err(super::super::Error::ConfigurationMissing {
                missing_config: "star_row_occurrence",
                ..
            })
        ));
        let row = Arc::new(Unused) as Arc<dyn StarRowAcquisition>;
        let output =
            session.batch_fixstars_acquired_ut_captured(&[None], f64::NAN, 0, false, &[Some(row)]);
        assert!(matches!(
            output.result,
            Err(super::super::Error::ConfigurationMissing {
                missing_config: "star_row_occurrence",
                ..
            })
        ));
    }
}
