//! House/time/sidereal composition over explicit session settings and images.
//!
//! The fixed-layout `houses*` methods return typed twelve-cusp records, or a
//! separate thirty-six-sector record for Gauquelin. The `public_houses*`
//! methods instead return tuples with variable-length cusp vectors and accept
//! `char` selectors. Neither family changes the session's configured observer;
//! geographic arguments belong to the current query.
use super::*;
use crate::houses::{
    HousesJdProviders, HousesSiderealProviders, fold_hsys_case, houses_ex_buffer_via,
    sidereal_selection_for_context,
};
use crate::sidereal::SealedAyanamshaServing;
use crate::time::{
    SiderealDeltaT, greenwich_apparent_sidereal_time_hours, house_armc_obliquity_deg,
};

#[cfg(test)]
mod live_house_anchor_tests {
    use super::*;
    use crate::stars::{
        AcquiredAnchorRow, AnchorRowAcquisition, AnchorSample, AnchorSampleOptions,
        NativeStarInputs,
    };
    use std::sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    };

    struct UnusedRows;
    impl AnchorRowAcquisition for UnusedRows {
        fn acquire_row(&self, _: &str) -> Result<Arc<dyn AcquiredAnchorRow>, Error> {
            panic!("unselected house operation read the live catalog")
        }
    }

    struct CachedRows(AtomicUsize);
    impl AnchorRowAcquisition for CachedRows {
        fn acquire_row(&self, key: &str) -> Result<Arc<dyn AcquiredAnchorRow>, Error> {
            assert_eq!(key, "SPICA");
            self.0.fetch_add(1, Ordering::SeqCst);
            Ok(Arc::new(CachedRow))
        }
    }
    struct CachedRow;
    impl AcquiredAnchorRow for CachedRow {
        fn sample(&self, _: f64, _: AnchorSampleOptions) -> Result<AnchorSample, Error> {
            Ok(AnchorSample::Cached(195.0))
        }
        fn galactic_inputs(&self) -> Result<NativeStarInputs, Error> {
            panic!("Spica must not acquire Galactic Center operands")
        }
    }

    fn session(mode: i32) -> Session {
        let mut ctx = EphemerisContext::new()
            .with_mode(ProviderMode::NativeJpl)
            .with_leb_file(None);
        ctx.set_sidereal_mode(mode, 0.0, 0.0);
        Session::from_context(ctx)
    }

    #[test]
    fn unconsumed_house_feeds_skip_catalog_main_and_site_validation() {
        let session = session(crate::constants::SIDM_TRUE_CITRA);
        let rows = UnusedRows;
        let result = session
            .with_native_houses_providers_with_position_context_with_anchor_rows(
                (999.0, 999.0, f64::NAN),
                &|_, _| panic!("unconsumed position resolver"),
                Some(&rows),
                |_, _, _, _| Ok(42),
            )
            .unwrap();
        assert_eq!(result, 42);
    }

    #[test]
    fn classical_house_offset_ignores_supplied_live_catalog() {
        let session = session(crate::constants::SIDM_FAGAN_BRADLEY);
        let rows = UnusedRows;
        session
            .with_native_houses_providers_with_anchor_rows(
                (0.0, 0.0, 0.0),
                Some(&rows),
                |_, _, sidereal, _| {
                    assert!((sidereal.mean_ut_plain_at)(2451545.0)?.is_finite());
                    assert!(
                        (sidereal.ayanamsha_ex_at)(2451545.0, FLG_NONUT)?
                            .1
                            .is_finite()
                    );
                    Ok(())
                },
            )
            .unwrap();
    }

    #[test]
    fn live_house_cache_is_consumed_only_by_selected_offset_callbacks() {
        let session = session(crate::constants::SIDM_TRUE_CITRA);
        let rows = CachedRows(AtomicUsize::new(0));
        session
            .with_native_houses_providers_with_anchor_rows(
                (0.0, 0.0, 0.0),
                Some(&rows),
                |_, _, sidereal, _| {
                    assert_eq!(rows.0.load(Ordering::SeqCst), 0);
                    assert_eq!((sidereal.mean_ut_plain_at)(2451545.0)?, 15.0);
                    assert_eq!((sidereal.ayanamsha_ex_at)(2451545.0, FLG_NONUT)?.1, 15.0);
                    assert_eq!(rows.0.load(Ordering::SeqCst), 2);
                    Ok(())
                },
            )
            .unwrap();
    }

    #[test]
    fn house_anchor_refusal_follows_original_clock_prefix() {
        struct Refused(AtomicUsize);
        impl AnchorRowAcquisition for Refused {
            fn acquire_row(&self, _: &str) -> Result<Arc<dyn AcquiredAnchorRow>, Error> {
                self.0.fetch_add(1, Ordering::SeqCst);
                Err(Error::ConfigurationMissing {
                    missing_config: "actual house anchor refusal",
                    suggestion: "keep the owned callback error",
                })
            }
        }
        let session = session(crate::constants::SIDM_TRUE_CITRA);
        let rows = Refused(AtomicUsize::new(0));
        session
            .with_native_houses_providers_with_anchor_rows(
                (0.0, 0.0, 0.0),
                Some(&rows),
                |_, _, sidereal, _| {
                    assert!((sidereal.mean_ut_plain_at)(f64::NAN).is_err());
                    assert_eq!(rows.0.load(Ordering::SeqCst), 0);
                    assert!(matches!(
                        (sidereal.mean_ut_plain_at)(2451545.0),
                        Err(Error::ConfigurationMissing {
                            missing_config: "actual house anchor refusal",
                            ..
                        })
                    ));
                    assert_eq!(rows.0.load(Ordering::SeqCst), 1);
                    Ok(())
                },
            )
            .unwrap();
    }
}

/// Synchronous configuration resolver borrowed by consumed house positions.
///
/// The shared context is immutable and the returned snapshot owns its registry.
/// The callable is borrowed only for the current operation; it retains no input
/// reference and supplies no astronomical position or clock value.
type HousePositionResolver<'a> =
    dyn Fn(&EphemerisContext, i32) -> Result<EphemerisContext, Error> + 'a;

/// Borrow one optional nutation acquisition and consume only missing words.
/// Inputs use UT1 Julian days and degree angles. A supplied NaN or signed zero
/// is not replaced. When both values exist the acquisition callback is unused;
/// otherwise its error precedes `run`. Callbacks are synchronous and no borrow
/// survives the call. This helper owns selection/order, not a time formula.
fn with_public_sidereal_components(
    jd_ut: f64,
    obliquity: Option<f64>,
    nutation: Option<f64>,
    acquire: &dyn Fn(f64) -> Result<[f64; 6], Error>,
    run: &dyn Fn(f64, f64) -> Result<f64, Error>,
) -> Result<f64, Error> {
    let mut eps = obliquity;
    let mut nut = nutation;
    if eps.is_none() || nut.is_none() {
        let words = acquire(jd_ut)?;
        if eps.is_none() {
            eps = Some(words[0]);
        }
        if nut.is_none() {
            nut = Some(words[2]);
        }
    }
    run(
        eps.expect("selected obliquity"),
        nut.expect("selected nutation"),
    )
}

impl Session {
    /// Borrow native house acquisition feeds for one synchronous operation.
    ///
    /// The session must explicitly select NativeJpl. Its owned registry is
    /// borrowed for DE-backed live sidereal anchors through the context star
    /// factory; construction performs no body sample. House frames use the
    /// configured Delta-T chain, while
    /// ayanamsha observation TT uses the independent bundled UT1 timescale.
    /// Extended UT ayanamsha retains the public UT→TT→UT entry conversion
    /// before the mean/true helper reacquires its scalar observation TT.
    /// USER_UT and Valens defining epochs retain the configured model and are
    /// acquired only by their selected branches. Sunshine consumes the genuine
    /// public equatorial Sun calc_ut route. Gauquelin body places similarly use
    /// public calc_ut; a per-call site replaces the observer only when TOPOCTR
    /// is consumed. Neither the caller's context nor its registry is mutated.
    ///
    /// `site` is `(longitude_deg_east, latitude_deg, height_m)`; its validity
    /// is checked only for a consumed topocentric body call. Gauquelin methods
    /// 2–5 use the genuine public rise/set search at each original UT1 query,
    /// with a zero-degree horizon, that query's observer and atmosphere, and
    /// disc/refraction bits selected by the sector kernel. Ordinary events
    /// validate a temporary observer even when the caller's flags omit TOPOCTR.
    /// The event search retains its own configured position/frame clocks and
    /// independent bundled sidereal clock; no epoch is inverted by this feed.
    /// A no-event result is preserved for the sector kernel to reject. A failed
    /// initial out-of-map positionability probe and bare math-domain errors map
    /// to GauquelinRiseSetNotFound; coverage, coordinate, configuration and
    /// later lookup errors preserve their typed identity. Methods 0/1 consume
    /// no event feed. All provider/image borrows end with `run`.
    /// This SDK form uses only the session's explicit configuration. Bindings
    /// can supply lazy position configuration with
    /// [`Self::with_native_houses_providers_with_position_context`].
    /// The existing native-ray callback owns catalog propagation and apparent
    /// reductions, including vector NOABERR, at each scalar TT node; the shared
    /// biased mean frame and eight-day interpolation are retained. Native anchor
    /// deflector coverage is strict at observation and closest
    /// approach. Earth/interpolation-node and deflector raw range failures
    /// become public ranges at the original anchor request TT, using the
    /// active main Sun-image envelope and absent star-body identity.
    pub fn with_native_houses_providers<T>(
        &self,
        site: (f64, f64, f64),
        run: impl FnOnce(
            crate::sidereal::SiderealSelection,
            &HousesJdProviders<'_>,
            &crate::houses::HousesSiderealCallbacks<'_>,
            &crate::houses::GauquelinSectorProviders<'_>,
        ) -> Result<T, Error>,
    ) -> Result<T, Error> {
        self.native_houses_providers(site, None, None, None, None, run)
    }

    /// Borrow native house feeds with call-owned live sidereal catalog rows.
    ///
    /// Feed signatures, units, frame/model clocks and errors follow
    /// [`Self::with_native_houses_providers`]. `rows` is consumed only when the
    /// selected house operation evaluates a dynamic ayanamsha; tropical,
    /// classical, fixed-epoch and unselected callbacks read no catalog. The
    /// released `houses_ex` offset occurs after tropical and projection routing.
    /// Every consumed helper/node uses the supplied row/cache owner, with lazy
    /// independent main acquisition on a miss. No row survives by substitution
    /// after refusal. `None` preserves the standalone static catalog behavior.
    pub fn with_native_houses_providers_with_anchor_rows<T>(
        &self,
        site: (f64, f64, f64),
        rows: Option<&dyn crate::stars::AnchorRowAcquisition>,
        run: impl FnOnce(
            crate::sidereal::SiderealSelection,
            &HousesJdProviders<'_>,
            &crate::houses::HousesSiderealCallbacks<'_>,
            &crate::houses::GauquelinSectorProviders<'_>,
        ) -> Result<T, Error>,
    ) -> Result<T, Error> {
        self.native_houses_providers(site, None, rows, None, None, run)
    }

    /// Borrow native house feeds with lazy public-position configuration.
    ///
    /// This is the configuration-boundary form of
    /// [`Self::with_native_houses_providers`]. `resolve` is called only when a
    /// Sunshine Sun, Gauquelin body or Gauquelin event callback consumes a
    /// genuine [`Self::calc_ut`] result. It receives the selected owned context and
    /// public body identifier; a consumed topocentric site is validated and
    /// installed before resolution (except the event positionability probe,
    /// which precedes observer validation). The returned context owns its registry
    /// snapshot. Its original UT1 Julian day and position flags are forwarded
    /// unchanged to the public position entry.
    ///
    /// Ordinary house frames and unselected Sun/body branches never call
    /// `resolve`. Live sidereal anchors continue to borrow this session's
    /// registry and strict native-ray provider. Configuration errors propagate
    /// before position acquisition; numerical and provider errors keep the
    /// same downstream house fallback rules. Resolver and provider borrows end
    /// with `run`; neither the caller's context nor its observer is modified.
    pub fn with_native_houses_providers_with_position_context<T>(
        &self,
        site: (f64, f64, f64),
        resolve: &dyn Fn(&EphemerisContext, i32) -> Result<EphemerisContext, Error>,
        run: impl FnOnce(
            crate::sidereal::SiderealSelection,
            &HousesJdProviders<'_>,
            &crate::houses::HousesSiderealCallbacks<'_>,
            &crate::houses::GauquelinSectorProviders<'_>,
        ) -> Result<T, Error>,
    ) -> Result<T, Error> {
        self.native_houses_providers(site, Some(resolve), None, None, None, run)
    }

    /// Borrow live house anchors alongside lazy public-position configuration.
    ///
    /// The position resolver has exactly the Sun/body/event roles of
    /// [`Self::with_native_houses_providers_with_position_context`]; it is not
    /// an anchor-row or numerical callback. `rows` has the lazy, synchronous
    /// catalog ownership of [`Self::with_native_houses_providers_with_anchor_rows`].
    /// Frames, site guards and both position and anchor clocks retain their
    /// separate stages. Neither borrowed capability escapes `run`.
    pub fn with_native_houses_providers_with_position_context_with_anchor_rows<T>(
        &self,
        site: (f64, f64, f64),
        resolve: &dyn Fn(&EphemerisContext, i32) -> Result<EphemerisContext, Error>,
        rows: Option<&dyn crate::stars::AnchorRowAcquisition>,
        run: impl FnOnce(
            crate::sidereal::SiderealSelection,
            &HousesJdProviders<'_>,
            &crate::houses::HousesSiderealCallbacks<'_>,
            &crate::houses::GauquelinSectorProviders<'_>,
        ) -> Result<T, Error>,
    ) -> Result<T, Error> {
        self.native_houses_providers(site, Some(resolve), rows, None, None, run)
    }

    /// Assemble borrowed house feeds without consuming unselected callbacks.
    ///
    /// `site` uses east-positive longitude/latitude degrees and height metres.
    /// The optional resolver belongs only to consumed public Sun/body positions
    /// and genuine event-search position calls. Event positionability probes
    /// precede observer validation; ordinary topocentric samples follow it. Frame/model and strict anchor feeds
    /// retain the caller's context. An optional model clock supplies only
    /// model TT and selected defining epochs; independent anchor time remains
    /// unchanged. The warning observer runs only at a consumed unknown-mode
    /// fallback. All closures execute synchronously inside `run`; their image,
    /// resolver and local-context borrows cannot escape.
    pub(super) fn native_houses_providers<T>(
        &self,
        site: (f64, f64, f64),
        resolve: Option<&HousePositionResolver<'_>>,
        rows: Option<&dyn crate::stars::AnchorRowAcquisition>,
        model_delta: Option<&dyn Fn(f64) -> Result<f64, Error>>,
        unknown_mode_at: Option<&dyn Fn(i32) -> Result<(), Error>>,
        run: impl FnOnce(
            crate::sidereal::SiderealSelection,
            &HousesJdProviders<'_>,
            &crate::houses::HousesSiderealCallbacks<'_>,
            &crate::houses::GauquelinSectorProviders<'_>,
        ) -> Result<T, Error>,
    ) -> Result<T, Error> {
        if self.ctx.mode() != ProviderMode::NativeJpl {
            return Err(Error::ConfigurationMissing {
                missing_config: "native_houses_mode",
                suggestion: "select NativeJpl for native house providers",
            });
        }
        crate::positions::with_context_star_serving(
            &self.ctx,
            self.delta_t_model(),
            self.ctx.effective_iers_delta_t(),
            0,
            None,
            |star_serving, _| {
                let anchors = star_serving.anchor_provider();
                let (selection, _) = sidereal_selection_for_context(&self.ctx);
                let base = crate::sidereal::strip_to_base_mode(selection.mode());
                let rows = rows.filter(|_| {
                    crate::sidereal::DYNAMIC_AYANAMSHA_MODES.contains(&base)
                        && !matches!(
                            base,
                            crate::constants::SIDM_GALALIGN_MARDYKS
                                | crate::constants::SIDM_VALENS_MOON
                        )
                });
                let delta = |ut| match model_delta {
                    Some(acquire) => acquire(ut),
                    None => self.deltat(ut),
                };
                let aya = SealedAyanamshaServing::new(selection, &delta, Some(&anchors))
                    .with_unknown_mode_observer(unknown_mode_at);
                let frame = |ut| match model_delta {
                    Some(acquire) => self.house_frame_with_delta(ut, acquire),
                    None => self.frame_at(ut),
                };
                let sun = |ut, eph| {
                    native_house_position_with_context(
                        self,
                        ut,
                        SUN,
                        FLG_EQUATORIAL | (eph & (FLG_JPLEPH | FLG_SWIEPH)),
                        resolve,
                    )
                    .map(|p| p.values[1])
                };
                let extended = |ut, masked| {
                    if let Some(rows) = rows {
                        return self.house_live_extended_offset(ut, masked, selection, rows);
                    }
                    let entry_tt = crate::time::default_timescale_tt_from_ut1(ut)?;
                    let helper_ut = crate::time::default_timescale_ut1_from_tt(entry_tt)?;
                    let evaluate = |aya: &SealedAyanamshaServing<'_>| {
                        crate::planetary_moons::moon_spk_ayanamsha(
                            helper_ut,
                            masked & (FLG_TRUEPOS | FLG_NOABERR | FLG_NOGDEFL | FLG_NONUT),
                            aya,
                        )
                    };
                    let value = evaluate(&aya)?;
                    let base = crate::sidereal::strip_to_base_mode(selection.mode());
                    Ok((crate::sidereal::ayanamsha_ex_retflag(masked, base), value))
                };
                let mean = |ut| {
                    if let Some(rows) = rows {
                        return self.house_live_mean_offset(ut, selection, rows);
                    }
                    let tt = crate::time::default_timescale_tt_from_ut1(ut)?;
                    aya.mean_tt_with(tt, false, false).map(|(value, _)| value)
                };
                let nutation = |tt| crate::erfa::nut06a(tt);
                let user = || aya.user_pair_tt().map(|(tt, _)| tt);
                let valens = || aya.valens_epoch_tt();
                let body = |ut, body_id, flags| {
                    let position = if flags & FLG_TOPOCTR != 0 {
                        let observer = Observer::new(site.0, site.1, site.2)?;
                        let selected = Self::from_context(self.ctx.clone().with_observer(observer));
                        native_house_position_with_context(&selected, ut, body_id, flags, resolve)?
                    } else {
                        native_house_position_with_context(self, ut, body_id, flags, resolve)?
                    };
                    Ok((position.values[0], position.values[1]))
                };
                let jd = HousesJdProviders {
                    frame_at: &frame,
                    sun_declination_at: &sun,
                    deltat_days_at: &delta,
                };
                let sidereal = crate::houses::HousesSiderealCallbacks {
                    ayanamsha_ex_at: &extended,
                    mean_ut_plain_at: &mean,
                    nutation_at: &nutation,
                    user_t0_tt_at: &user,
                    valens_epoch_tt_at: &valens,
                };
                let rise_set = |query| native_gauquelin_rise_set(self, query, resolve);
                let gauquelin = crate::houses::GauquelinSectorProviders {
                    frame_at: &frame,
                    body_ecliptic_at: &body,
                    rise_set_at: Some(&rise_set),
                };
                run(selection, &jd, &sidereal, &gauquelin)
            },
        )
    }

    /// Live house offset through the released extended UT clock prefix.
    /// Original UT1 acquires independent TT then helper UT1 before any row,
    /// cache or image operation. The helper reacquires observation TT inside
    /// the accepted mean/true kernel; configured Delta-T remains available
    /// only for selected defining epochs. Caller flags retain their old mask.
    pub(super) fn house_live_extended_offset(
        &self,
        ut: f64,
        masked: i32,
        selection: crate::sidereal::SiderealSelection,
        rows: &dyn crate::stars::AnchorRowAcquisition,
    ) -> Result<(i32, f64), Error> {
        let entry_tt = crate::time::default_timescale_tt_from_ut1(ut)?;
        let helper_ut = crate::time::default_timescale_ut1_from_tt(entry_tt)?;
        let delta = |ut| self.deltat(ut);
        let value = self
            .star_anchor_provider_captured(rows, |anchors| {
                let aya = SealedAyanamshaServing::new(selection, &delta, Some(anchors));
                crate::planetary_moons::moon_spk_ayanamsha(
                    helper_ut,
                    masked & (FLG_TRUEPOS | FLG_NOABERR | FLG_NOGDEFL | FLG_NONUT),
                    &aya,
                )
            })
            .result?;
        let base = crate::sidereal::strip_to_base_mode(selection.mode());
        Ok((crate::sidereal::ayanamsha_ex_retflag(masked, base), value))
    }

    /// Live geometric projection zero point at independent observation TT.
    /// The original UT1-to-TT prefix precedes catalog/cache/main acquisition.
    /// Every helper/node retains the same call-owned supplier; no fast-reader
    /// mean entry is substituted for this defining projection calculation.
    pub(super) fn house_live_mean_offset(
        &self,
        ut: f64,
        selection: crate::sidereal::SiderealSelection,
        rows: &dyn crate::stars::AnchorRowAcquisition,
    ) -> Result<f64, Error> {
        let tt = crate::time::default_timescale_tt_from_ut1(ut)?;
        let delta = |ut| self.deltat(ut);
        self.star_anchor_provider_captured(rows, |anchors| {
            SealedAyanamshaServing::new(selection, &delta, Some(anchors))
                .mean_tt_with(tt, false, false)
                .map(|(value, _)| value)
        })
        .result
    }

    /// Checked canonical clock: consume only the requested TT-minus-UT1.
    /// The legacy slots of this private carrier are inert, so no unrelated
    /// boundary clock is acquired. Actual public legacy offsets/getters and
    /// serving-constructor contracts remain with their explicit utilities.
    pub(super) fn sidereal_delta(&self, jd_ut: f64) -> Result<SiderealDeltaT, Error> {
        SiderealDeltaT::try_new(self.deltat(jd_ut)?, 0.0, 0.0)
    }

    /// Producer extension: `(ARMC at longitude zero, true obliquity)` in
    /// degrees at UT, using the configured body/model clock and native nut06a.
    /// Used by houses and caller geometry; does not evaluate an unbound body.
    pub fn frame_at(&self, jd_ut: f64) -> Result<(f64, f64), Error> {
        let dt = self.sidereal_delta(jd_ut)?;
        let jd_tt = jd_ut + self.deltat(jd_ut)?;
        let (dpsi, deps) = crate::erfa::nut06a(jd_tt)?;
        house_armc_obliquity_deg(jd_ut, &dt, dpsi.to_degrees(), deps.to_degrees())
    }

    /// House frame from the actual caller-model clock, in its consumption
    /// order. Original UT1 first acquires TT for nutation/obliquity, then the
    /// independently consumed sidereal TT. Both consume the requested epoch;
    /// no boundary calibration enters canonical time. Angles are degrees;
    /// each clock/provider error stops before the next acquisition.
    pub(super) fn house_frame_with_delta(
        &self,
        jd_ut: f64,
        delta: &dyn Fn(f64) -> Result<f64, Error>,
    ) -> Result<(f64, f64), Error> {
        let obliquity_tt = jd_ut + delta(jd_ut)?;
        let (dpsi, deps) = crate::erfa::nut06a(obliquity_tt)?;
        let eps_true =
            crate::positions::vondrak_mean_obliquity_deg(obliquity_tt)? + deps.to_degrees();
        let gmst = crate::time::public_sidereal::mean_deg(jd_ut, delta)?;
        let armc = (gmst + dpsi.to_degrees() * eps_true.to_radians().cos()).rem_euclid(360.);
        Ok((armc, eps_true))
    }

    /// Greenwich apparent sidereal hours at UT1 Julian Day.
    /// Native mode first consumes public obliquity/nutation, then the selected
    /// sidereal model clocks. Stored modes retain their checked house-frame
    /// contract. Values use a radians/degrees round-trip and modulo24.
    pub fn sidtime(&self, jd_ut: f64) -> Result<f64, Error> {
        if self.ctx.mode() == ProviderMode::NativeJpl {
            let nutation = self.calc_ut(jd_ut, ECL_NUT, 0)?;
            return self.sidtime0(jd_ut, nutation.values[0], nutation.values[2]);
        }
        let dt = self.sidereal_delta(jd_ut)?;
        let (_, eps) = self.frame_at(jd_ut)?;
        let (dpsi, _) = crate::erfa::nut06a(jd_ut + self.deltat(jd_ut)?)?;
        greenwich_apparent_sidereal_time_hours(jd_ut, eps, dpsi.to_degrees(), &dt)
    }

    /// Local apparent sidereal hours at UT1 with optional true obliquity and
    /// nutation in longitude, both in degrees. Longitude is east-positive.
    /// If either component is absent, one public calc_ut(ECL_NUT,0) supplies
    /// only the missing component; a supplied value keeps its original word.
    /// This acquisition precedes the selected sidereal model clock. Longitude
    /// is added before modulo24. Native raw IEEE and stored checked domains
    /// remain distinct; no image, model or extra sample is acquired eagerly.
    pub fn sidtime_with(
        &self,
        jd_ut: f64,
        longitude: f64,
        obliquity: Option<f64>,
        nutation: Option<f64>,
    ) -> Result<f64, Error> {
        with_public_sidereal_components(
            jd_ut,
            obliquity,
            nutation,
            &|ut| self.calc_ut(ut, ECL_NUT, 0).map(|position| position.values),
            &|eps, nut| {
                if self.ctx.mode() == ProviderMode::NativeJpl {
                    return crate::time::public_sidereal::hours_with_longitude(
                        jd_ut,
                        longitude,
                        eps,
                        nut,
                        &|ut| self.deltat(ut),
                    );
                }
                crate::time::local_apparent_sidereal_time_hours(
                    jd_ut,
                    longitude,
                    eps,
                    nut,
                    &self.sidereal_delta(jd_ut)?,
                )
            },
        )
    }

    /// Sidereal hours at UT1 with supplied true obliquity/nutation in degrees.
    /// Native mode consumes only the selected branch's model Delta-T calls:
    /// only the requested epoch, without boundary calibration. NaN model
    /// offsets and angles propagate; infinite cosine
    /// arguments fail after the clock. Stored modes retain finite checks.
    pub fn sidtime0(&self, jd_ut: f64, eps: f64, nut: f64) -> Result<f64, Error> {
        if self.ctx.mode() == ProviderMode::NativeJpl {
            return crate::time::public_sidereal::hours(jd_ut, eps, nut, &|ut| self.deltat(ut));
        }
        greenwich_apparent_sidereal_time_hours(jd_ut, eps, nut, &self.sidereal_delta(jd_ut)?)
    }

    /// Mean ayanamsha in degrees at TT, retaining the established TT inversion.
    pub fn get_ayanamsa(&self, jd_tt: f64) -> Result<f64, Error> {
        let (selection, _) = sidereal_selection_for_context(&self.ctx);
        let delta = |ut| self.deltat(ut);
        SealedAyanamshaServing::new(selection, &delta, None)
            .mean_tt(jd_tt)
            .map(|(v, _)| v)
    }

    /// Extended TT ayanamsha: degrees, returned flags and producer advisory.
    /// Unsupported dynamic anchors refuse instead of substituting a constant.
    pub fn get_ayanamsa_ex(&self, jd_tt: f64, flags: i32) -> Result<Ayanamsha, Error> {
        let (selection, _) = sidereal_selection_for_context(&self.ctx);
        let delta = |ut| self.deltat(ut);
        let (returned_flags, value, warning) =
            SealedAyanamshaServing::new(selection, &delta, None).ex_tt(jd_tt, flags)?;
        Ok(Ayanamsha {
            value,
            returned_flags,
            diagnostic: warning.map_or_else(String::new, |w| w.to_string()),
        })
    }

    /// Mean ayanamsha in degrees at UT, preserving selected model/USER epochs.
    pub fn get_ayanamsa_ut(&self, jd_ut: f64) -> Result<f64, Error> {
        let (selection, _) = sidereal_selection_for_context(&self.ctx);
        let delta = |jd| self.deltat(jd);
        SealedAyanamshaServing::new(selection, &delta, None)
            .mean_ut(jd_ut)
            .map(|(v, _)| v)
    }

    /// Extended flag-aware ayanamsha in degrees at UT, with returned flags
    /// and any producer advisory. Missing required anchors remain errors.
    pub fn get_ayanamsa_ex_ut(&self, jd_ut: f64, flags: i32) -> Result<Ayanamsha, Error> {
        let (selection, _) = sidereal_selection_for_context(&self.ctx);
        let delta = |jd| self.deltat(jd);
        let (returned_flags, value, warning) =
            SealedAyanamshaServing::new(selection, &delta, None).ex_ut(jd_ut, flags)?;
        Ok(Ayanamsha {
            value,
            returned_flags,
            diagnostic: warning.map_or_else(String::new, |w| w.to_string()),
        })
    }

    /// Twelve tropical cusps and eight angles at a UT1 Julian day, in degrees.
    ///
    /// `lat` is north-positive geographic latitude in `[-90, 90]` degrees;
    /// `lon` is east-positive longitude in `[-180, 360]` degrees. `hsys` is an
    /// ASCII byte such as `b'P'` (Placidus), `b'O'` (Porphyry) or `b'W'`
    /// (whole sign); see [`crate::houses`] for the system inventory.
    ///
    /// Returns [`Houses`] with house one at index zero and no padding. This
    /// is [`Self::houses_ex`] with flags zero, including its fixed-layout
    /// refusals and errors. No alternate system is selected on a polar error;
    /// use [`Self::public_houses_with_fallback`] for an explicit fallback.
    pub fn houses(&self, jd_ut: f64, lat: f64, lon: f64, hsys: u8) -> Result<Houses, Error> {
        self.houses_ex(jd_ut, 0, lat, lon, hsys)
    }

    /// Calculate a fixed twelve-cusp layout with explicit frame/output flags.
    ///
    /// # Arguments
    ///
    /// - `jd_ut`: UT1 Julian day, not UTC or TT.
    /// - `flags`: `FLG_SIDEREAL` selects the context's ayanamsha/projection;
    ///   `FLG_RADIANS` converts central positions to radians. Ephemeris bits
    ///   affect the Sunshine Sun acquisition, not selection of a different
    ///   context provider. Other body-position flags are not a general house
    ///   coordinate interface.
    /// - `lat`, `lon`: finite north/east degrees in `[-90, 90]` and
    ///   `[-180, 360]`, respectively.
    /// - `hsys`: ASCII alphabetic house-system byte; known codes are listed in
    ///   [`crate::houses`]. Ordinary lowercase codes fold to uppercase, while
    ///   `i` remains the alternate Sunshine construction. Unknown letters use
    ///   the kernel's Placidus-default construction.
    ///
    /// # Returns
    ///
    /// [`Houses`]: twelve cusps and eight angles, in degrees unless RADIANS is
    /// requested. The default is tropical geometry of date; sidereal flags
    /// select the configured zero point or projection plane. This method
    /// rejects both `G` and `g`: use [`Self::houses_gauquelin`] for a fixed
    /// thirty-six-sector record or [`Self::public_houses_ex`] for flexible
    /// `char`/vector layouts. Geographic inputs do not mutate the observer.
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnsupportedFlagCombination`] for Gauquelin or a
    /// nonalphabetic byte, [`Error::InvalidObserver`] for invalid coordinates,
    /// and clock/frame/numerical errors for invalid consumed epochs/geometry.
    /// [`Error::PolarCircle`] propagates when the selected construction is
    /// undefined at the site; no alternate system is automatically selected.
    /// Required Sunshine or sidereal-anchor data/configuration errors retain
    /// their provider categories. Most tropical systems need no planetary file.
    pub fn houses_ex(
        &self,
        jd_ut: f64,
        flags: i32,
        lat: f64,
        lon: f64,
        hsys: u8,
    ) -> Result<Houses, Error> {
        self.house_sample(jd_ut, flags, lat, lon, hsys)
            .map(|(cusps, angles)| Houses { cusps, angles })
    }

    /// Return fixed-layout house positions and rates per UT1 day.
    ///
    /// Arguments, frames, selector restrictions and errors follow
    /// [`Self::houses_ex`]. [`HousesWithSpeeds`] contains twelve cusp positions,
    /// eight angles, twelve cusp rates and eight angle rates, with matching
    /// zero-based indices. The centered/Richardson stencil samples the same
    /// house function at neighbouring UT1 epochs; side-sample errors propagate.
    ///
    /// `FLG_RADIANS` converts **central positions only**. Rates always remain
    /// degrees per UT1 day, not radians/day or rates per TT day. Existing wrap,
    /// sign-boundary and discontinuity conventions can return zero rates;
    /// a returned rate is not proof of smoothness at a branch boundary.
    /// Both `G` and `g` are refused before acquisition, not truncated.
    pub fn houses_ex2(
        &self,
        jd_ut: f64,
        flags: i32,
        lat: f64,
        lon: f64,
        hsys: u8,
    ) -> Result<HousesWithSpeeds, Error> {
        if hsys.eq_ignore_ascii_case(&b'G') {
            return Err(Error::UnsupportedFlagCombination {
                entry: "api::houses_ex2 (twelve cusps)",
                flags,
            });
        }
        // The existing stencil wraps degrees, not radians. Shape only the
        // central positions after differentiation; rates retain degrees/day.
        let sample = |jd| self.house_sample_buffer(jd, flags & !FLG_RADIANS, lat, lon, hsys);
        let (mut cusps, mut angles, cusp_speeds, angle_speeds) =
            crate::houses::houses_ex2_buffer_via(&sample, jd_ut, fold_hsys_case(char::from(hsys)))?;
        if flags & FLG_RADIANS != 0 {
            cusps
                .iter_mut()
                .for_each(|value| *value = value.to_radians());
            angles = angles.map(f64::to_radians);
        }
        Ok(HousesWithSpeeds {
            cusps: twelve(cusps)?,
            angles,
            cusp_speeds: twelve(cusp_speeds)?,
            angle_speeds,
        })
    }

    /// Return all thirty-six Gauquelin sectors and their rates at UT1.
    ///
    /// `jd_ut`, `flags`, `lat` and `lon` have the units and frame meanings of
    /// [`Self::houses_ex`]. The result is [`GauquelinSectors`], with sector one
    /// at index zero, eight angles in [`Houses::angles`] order, and matching
    /// rates. RADIANS converts central positions only; all rates are degrees
    /// per UT1 day. Errors, side-sample domains and discontinuity limitations
    /// follow [`Self::houses_ex2`]. Polar geometry errors propagate; this method
    /// neither truncates sectors nor substitutes another house system.
    pub fn houses_gauquelin(
        &self,
        jd_ut: f64,
        flags: i32,
        lat: f64,
        lon: f64,
    ) -> Result<GauquelinSectors, Error> {
        let sample = |jd| self.house_sample_buffer(jd, flags & !FLG_RADIANS, lat, lon, b'G');
        let (mut sectors, mut angles, speeds, angle_speeds) =
            crate::houses::houses_ex2_buffer_via(&sample, jd_ut, 'G')?;
        if flags & FLG_RADIANS != 0 {
            sectors
                .iter_mut()
                .for_each(|value| *value = value.to_radians());
            angles = angles.map(f64::to_radians);
        }
        let convert = |v: crate::houses::CuspBuffer| {
            v.as_slice()
                .try_into()
                .map_err(|_| Error::ConfigurationMissing {
                    missing_config: "gauquelin_layout",
                    suggestion: "expected thirty-six sectors",
                })
        };
        Ok(GauquelinSectors {
            sectors: convert(sectors)?,
            angles,
            sector_speeds: convert(speeds)?,
            angle_speeds,
        })
    }

    /// Twelve-cusp projection with layout validation, never unchecked indexing.
    fn house_sample(
        &self,
        jd: f64,
        flags: i32,
        lat: f64,
        lon: f64,
        hsys: u8,
    ) -> Result<([f64; 12], [f64; 8]), Error> {
        if hsys.eq_ignore_ascii_case(&b'G') {
            return Err(Error::UnsupportedFlagCombination {
                entry: "api::houses (twelve cusps)",
                flags,
            });
        }
        let (cusps, angles) = self.house_sample_buffer(jd, flags, lat, lon, hsys)?;
        Ok((twelve(cusps)?, angles))
    }

    /// Assemble one houses sample; all provider borrows end in this call.
    /// NativeJpl uses the selected native callbacks: configured UT-frame TT,
    /// independent anchor observation TT, and genuine public Sun acquisition.
    /// Other modes retain their prior eager model/provider assembly.
    fn house_sample_buffer(
        &self,
        jd: f64,
        flags: i32,
        lat: f64,
        lon: f64,
        hsys: u8,
    ) -> Result<(crate::houses::CuspBuffer, [f64; 8]), Error> {
        if !hsys.is_ascii_alphabetic() {
            return Err(Error::UnsupportedFlagCombination {
                entry: "api::houses (twelve cusps)",
                flags,
            });
        }
        if self.ctx.mode() == ProviderMode::NativeJpl {
            return self.with_native_houses_providers(
                (0.0, 0.0, 0.0),
                |selection, jd_providers, sid_providers, _| {
                    crate::houses::houses_ex_buffer_via_callbacks(
                        jd,
                        lat,
                        lon,
                        char::from(hsys),
                        flags,
                        selection,
                        jd_providers,
                        sid_providers,
                    )
                },
            );
        }
        let delta = |jd| self.deltat(jd);
        let frame = |jd| self.frame_at(jd);
        let sun = |jd, eph| {
            self.calc_ut(
                jd,
                SUN,
                (eph & (FLG_JPLEPH | FLG_SWIEPH)) | FLG_EQUATORIAL | (flags & FLG_TOPOCTR),
            )
            .map(|p| p.values[1])
        };
        let jd_providers = HousesJdProviders {
            frame_at: &frame,
            sun_declination_at: &sun,
            deltat_days_at: &delta,
        };
        let (selection, _) = sidereal_selection_for_context(&self.ctx);
        let ayanamsha = SealedAyanamshaServing::new(selection, &delta, None);
        let selection = ayanamsha.selection();
        let extended = |jd, masked| {
            ayanamsha
                .ex_ut(jd, masked)
                .map(|(ret, value, _)| (ret, value))
        };
        let mean = |jd| ayanamsha.mean_ut(jd).map(|(value, _)| value);
        let nutation = |tt| crate::erfa::nut06a(tt);
        let (_, t0, _) = selection.full();
        let user_t0_tt = if selection.bits() & SIDBIT_USER_UT != 0 {
            t0 + delta(t0)?
        } else {
            t0
        };
        let valens_epoch_tt = ayanamsha.valens_epoch_tt()?;
        let sid_providers = HousesSiderealProviders {
            ayanamsha_ex_at: &extended,
            mean_ut_plain_at: &mean,
            nutation_at: &nutation,
            user_t0_tt,
            valens_epoch_tt,
        };
        houses_ex_buffer_via(
            jd,
            lat,
            lon,
            char::from(hsys),
            flags,
            selection,
            &jd_providers,
            &sid_providers,
        )
    }
}

/// Test-only retained range-policy control for the original house draft.
/// Promote a raw native anchor range miss at its original scalar TT request.
/// The main Sun image supplies the full descriptor envelope, independently
/// of which Earth node or deflector/closest-approach read failed. Bounds are
/// Julian days from the existing image-metadata conversion; body is absent
/// for a star anchor. The raw error reaches this helper intact, and only
/// SpkOutOfRange is promoted. Other typed errors retain their identity.
/// This numeric contract does not reconstruct filenames/calendar messages.
#[cfg(test)]
fn native_anchor_range_error(
    provider: &crate::positions::LocalSpkStateProvider<'_>,
    request_tt: f64,
    error: Error,
) -> Error {
    if !matches!(&error, Error::SpkOutOfRange { .. }) {
        return error;
    }
    match provider.public_range_error(SUN, request_tt, error) {
        Error::EphemerisRange {
            requested_jd,
            start_jd,
            end_jd,
            ..
        } => Error::EphemerisRange {
            requested_jd,
            start_jd,
            end_jd,
            body: None,
        },
        other => other,
    }
}

/// Evaluate one public position after optional boundary configuration.
///
/// `ut` is the original UT1 Julian day and `flags` remains the caller's public
/// word. The resolver receives an immutable context, including any validated
/// temporary observer, and returns an owned snapshot whose images remain live
/// through `calc_ut`. No epoch inversion or body-state callback is introduced.
/// `None` preserves direct SDK configuration; resolver failures precede the
/// public position call and propagate unchanged.
fn native_house_position_with_context(
    session: &Session,
    ut: f64,
    body: i32,
    flags: i32,
    resolve: Option<&HousePositionResolver<'_>>,
) -> Result<Position, Error> {
    match resolve {
        None => session.calc_ut(ut, body, flags),
        Some(resolve) => {
            Session::from_context(resolve(&session.ctx, body)?).calc_ut(ut, body, flags)
        }
    }
}

/// Acquire one genuine rise/set result for an integer-body sector query.
///
/// The query carries an original UT1 Julian day, east-positive longitude and
/// north-positive latitude in degrees, height in metres, pressure in mbar,
/// temperature in Celsius and unchanged position flags. The public event
/// search receives those values exactly and a zero-degree horizon. Its
/// temporary observer, configured position/frame clocks, independent bundled
/// GAST clock, split provider reads and event solver belong to that search.
/// Return its status and first UT1 event component without shaping or retry;
/// `-2` and positive-zero event time remain the kernel's no-event input.
///
/// A synchronous configuration wrapper counts consumed public position calls,
/// never supplying coordinates or making a second probe. An out-of-map target
/// must first pass a positionability probe; its sole failed first lookup is
/// the illegal-direction refusal. A later ordinary lookup is not that refusal.
/// With no caller resolver the wrapper clones the unchanged immutable context
/// only when consumed. All context, registry and callback borrows end here;
/// neither the caller's observer nor its settings are modified.
///
/// Illegal-direction and bare math-domain failures are converted by
/// [`native_gauquelin_event_error`]. Every other typed failure propagates; no
/// geometric sector or alternate provider substitutes for a missing event.
fn native_gauquelin_rise_set(
    session: &Session,
    query: crate::houses::GauquelinRiseSetQuery,
    resolve: Option<&HousePositionResolver<'_>>,
) -> Result<(i32, f64), Error> {
    let position_calls = std::cell::Cell::new(0_usize);
    let selected_context = |context: &EphemerisContext, body| {
        position_calls.set(position_calls.get().saturating_add(1));
        match resolve {
            Some(resolve) => resolve(context, body),
            None => Ok(context.clone()),
        }
    };
    session
        .rise_trans_true_hor_with_position_context(
            query.jd_start_ut,
            query.body_id,
            query.rsmi,
            [query.lon_deg_east, query.lat_deg, query.height_m],
            query.atpress_mbar,
            query.attemp_c,
            0.0,
            query.flags,
            &selected_context,
        )
        .map(|(res, tret)| (res, tret[0]))
        .map_err(|error| native_gauquelin_event_error(&query, position_calls.get(), error))
}

/// Apply the sector's narrowly caught event-error policy without erasing ranges.
///
/// Disc-center/refraction bits identify the owning method 2–5. A math-domain
/// error denotes a bare value-domain failure. A matching UnknownBody from the
/// sole first consumed position call for an out-of-map target denotes the
/// event search's illegal-direction probe failure; that direction refusal is
/// caught here, while a later ordinary UnknownBody is not. Classical bodies
/// have no such probe. The event implementation must keep this phase contract.
///
/// Existing coverage errors retain their exact body, requested day and bounds.
/// Coordinate/undefined-direction, configuration, native-provider and explicit
/// nonfinite-input errors also pass through unchanged. Fixed-size geographic
/// and event records cannot produce sequence-length/index/coercion errors.
/// This policy uses typed variants and actual acquisition phase, never broad
/// error categories, display text or manufactured astronomical values.
fn native_gauquelin_event_error(
    query: &crate::houses::GauquelinRiseSetQuery,
    position_calls: usize,
    error: Error,
) -> Error {
    let classical = (SUN..=PLUTO).contains(&query.body_id) || query.body_id == EARTH;
    let caught = matches!(error.primary_error(), Error::MathDomain { .. })
        || (!classical
            && position_calls == 1
            && matches!(error.primary_error(), Error::UnknownBody { body_id } if *body_id == query.body_id));
    if !caught {
        return error;
    }
    let center = query.rsmi & BIT_DISC_CENTER != 0;
    let refracted = query.rsmi & BIT_NO_REFRACTION == 0;
    let method = match (center, refracted) {
        (true, false) => 2,
        (true, true) => 3,
        (false, false) => 4,
        (false, true) => 5,
    };
    Error::GauquelinRiseSetNotFound { method }.with_cause(error)
}

/// Convert a validated twelve-house output, refusing any unexpected layout.
fn twelve(values: impl AsRef<[f64]>) -> Result<[f64; 12], Error> {
    values
        .as_ref()
        .try_into()
        .map_err(|_| Error::ConfigurationMissing {
            missing_config: "house_layout",
            suggestion: "expected twelve house cusps",
        })
}

#[cfg(test)]
#[path = "native_houses_tests.rs"]
mod native_houses_tests;

#[cfg(test)]
mod public_time_binding_tests {
    use super::*;
    use std::cell::RefCell;

    #[test]
    fn supplied_components_use_no_nutation_feed() {
        let result = with_public_sidereal_components(
            2451545.0,
            Some(23.4),
            Some(-0.0),
            &|_| panic!("unused nutation feed"),
            &|eps, nut| {
                assert_eq!(eps.to_bits(), 23.4_f64.to_bits());
                assert_eq!(nut.to_bits(), (-0.0_f64).to_bits());
                Ok(7.0)
            },
        )
        .unwrap();
        assert_eq!(result, 7.0);
    }

    #[test]
    fn missing_obliquity_preserves_supplied_nutation_word() {
        let calls = RefCell::new(Vec::new());
        with_public_sidereal_components(
            2451545.0,
            None,
            Some(-0.0),
            &|jd| {
                calls.borrow_mut().push("nutation");
                assert_eq!(jd.to_bits(), 2451545.0_f64.to_bits());
                Ok([23.4, 0.0, 0.125, 0.0, 0.0, 0.0])
            },
            &|eps, nut| {
                calls.borrow_mut().push("hours");
                assert_eq!(eps.to_bits(), 23.4_f64.to_bits());
                assert_eq!(nut.to_bits(), (-0.0_f64).to_bits());
                Ok(7.0)
            },
        )
        .unwrap();
        assert_eq!(*calls.borrow(), vec!["nutation", "hours"]);
    }

    #[test]
    fn missing_nutation_preserves_supplied_nan_word() {
        let nan = f64::from_bits(0x7ff8000000000023);
        with_public_sidereal_components(
            2451545.0,
            Some(nan),
            None,
            &|_| Ok([23.4, 0.0, 0.125, 0.0, 0.0, 0.0]),
            &|eps, nut| {
                assert_eq!(eps.to_bits(), nan.to_bits());
                assert_eq!(nut.to_bits(), 0.125_f64.to_bits());
                Ok(7.0)
            },
        )
        .unwrap();
    }

    #[test]
    fn both_missing_share_one_actual_acquisition() {
        let calls = RefCell::new(0);
        with_public_sidereal_components(
            2451545.0,
            None,
            None,
            &|_| {
                *calls.borrow_mut() += 1;
                Ok([23.4, 0.0, 0.125, 0.0, 0.0, 0.0])
            },
            &|eps, nut| {
                assert_eq!(eps.to_bits(), 23.4_f64.to_bits());
                assert_eq!(nut.to_bits(), 0.125_f64.to_bits());
                Ok(7.0)
            },
        )
        .unwrap();
        assert_eq!(*calls.borrow(), 1);
    }

    #[test]
    fn acquisition_failure_precedes_hours_callback() {
        let error = with_public_sidereal_components(
            2451545.0,
            None,
            None,
            &|_| {
                Err(Error::ConfigurationMissing {
                    missing_config: "actual-test-feed",
                    suggestion: "test feed",
                })
            },
            &|_, _| panic!("hours reached after feed refusal"),
        )
        .unwrap_err();
        assert!(matches!(
            error,
            Error::ConfigurationMissing {
                missing_config: "actual-test-feed",
                ..
            }
        ));
    }
}
