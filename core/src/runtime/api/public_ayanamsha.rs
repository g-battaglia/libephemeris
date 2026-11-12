//! Scalar ayanamsha entries over explicit TT and UT1 observation clocks.
//!
//! One immutable session retains actual reader/main images. Reader-fast mean
//! UT1 uses stored Delta-T days; analytic observations use the bundled global
//! timescale. A TT entry retains its supplied epoch directly. Configured model Delta-T is consumed only by the delegated fast
//! inversion and USER_UT/Valens defining epochs. Angles are degrees.

use super::{Ayanamsha, Session};
use crate::constants::{
    FLG_NOABERR, FLG_NOGDEFL, FLG_NONUT, FLG_TRUEPOS, SIDBIT_ECL_DATE, SIDBIT_USER_UT,
    SIDM_GALALIGN_MARDYKS, SIDM_USER, SIDM_VALENS_MOON,
};
use crate::errors::{Error, Result, Warning};
use crate::sidereal::{
    DYNAMIC_AYANAMSHA_MODES, SealedAyanamshaServing, SiderealSelection, ayanamsha_ex_retflag,
    strip_to_base_mode,
};
use crate::stars::AnchorRowAcquisition;

#[path = "controlled_ayanamsha.rs"]
mod controlled_ayanamsha;

/// Catalog-free modes never acquire a main image or a star-serving stack.
/// Mardyks and Valens occupy the source dynamic set but use formula/epoch arms.
fn needs_anchor(selection: SiderealSelection) -> bool {
    let mode = strip_to_base_mode(selection.mode());
    DYNAMIC_AYANAMSHA_MODES.contains(&mode)
        && mode != SIDM_GALALIGN_MARDYKS
        && mode != SIDM_VALENS_MOON
}

/// Classify reader misses eligible for analytical sidereal evaluation.
/// Provisioning failures are handled by the caller. Classification uses typed
/// native errors, never message text.
fn fast_miss(error: &Error) -> bool {
    matches!(
        error.primary_error(),
        Error::LebOutOfRange { .. }
            | Error::LebBodyNotFound { .. }
            | Error::LebStarNotFound { .. }
            | Error::LebCorrupt { .. }
            | Error::MathDomain { .. }
    )
}

impl Session {
    /// Mean sidereal offset in degrees at an already acquired observation TT.
    ///
    /// `tt` is the scalar Terrestrial Time Julian day supplied by the actual
    /// observation-Time owner. `selection` retains the caller's mode, USER
    /// anchor and separately acquired projection bits. No reader Delta-T table,
    /// UT1 conversion or observation factory is consumed here. Only a selected
    /// USER_UT/Valens defining epoch consumes `configured_delta`, in days. Live
    /// anchor rows remain lazy and borrowed through this synchronous call.
    /// Errors and unknown-mode warnings retain the existing mean kernel's
    /// ownership and domain; this entry adds no alternate geometry or fallback.
    pub fn public_mean_at_observation_tt_with_selection(
        &self,
        tt: f64,
        selection: SiderealSelection,
        rows: Option<&dyn AnchorRowAcquisition>,
        configured_delta: &dyn Fn(f64) -> Result<f64>,
    ) -> Result<(f64, Option<Warning>)> {
        self.public_ayanamsha_selected_at_tt_with_delta(
            tt,
            0,
            false,
            selection,
            rows,
            configured_delta,
        )
    }

    /// Evaluate the source mean/true helper at one independent observation TT.
    /// Only a selected live catalog mode creates star serving. Its native main
    /// remains lazy until actual anchor acquisition; classical and unknown
    /// modes work without an astronomical image. The bound configured Delta-T
    /// closure is consumed only by USER_UT or Valens defining epochs.
    fn public_ayanamsha_at_tt(
        &self,
        tt: f64,
        flags: i32,
        extended: bool,
        rows: Option<&dyn AnchorRowAcquisition>,
    ) -> Result<(f64, Option<Warning>)> {
        let (selection, _) = crate::houses::sidereal_selection_for_context(&self.ctx);
        self.public_ayanamsha_selected_at_tt(tt, flags, extended, selection, rows)
    }

    /// Render a retained outer mode with the selected operation's current
    /// projection bits and USER epoch pair. The session owns the configuration
    /// and resources consumed after any caller-side flag operators.
    fn public_ayanamsha_selected_at_tt(
        &self,
        tt: f64,
        flags: i32,
        extended: bool,
        selection: SiderealSelection,
        rows: Option<&dyn AnchorRowAcquisition>,
    ) -> Result<(f64, Option<Warning>)> {
        let delta = |ut| self.deltat(ut);
        self.public_ayanamsha_selected_at_tt_with_delta(
            tt, flags, extended, selection, rows, &delta,
        )
    }

    /// Shared geometry with a call-owned configured defining-epoch clock.
    /// The supplier stays lazy: catalog-free modes without USER_UT/Valens do
    /// not invoke it, and no borrow escapes the synchronous mean reduction.
    fn public_ayanamsha_selected_at_tt_with_delta(
        &self,
        tt: f64,
        flags: i32,
        extended: bool,
        selection: SiderealSelection,
        rows: Option<&dyn AnchorRowAcquisition>,
        delta: &dyn Fn(f64) -> Result<f64>,
    ) -> Result<(f64, Option<Warning>)> {
        let noaberr = extended && flags & (FLG_TRUEPOS | FLG_NOABERR) != 0;
        let nogdefl = extended && flags & (FLG_TRUEPOS | FLG_NOGDEFL) != 0;
        let render = |aya: &SealedAyanamshaServing<'_>| {
            let (mean, warning) = aya.mean_tt_with(tt, noaberr, nogdefl)?;
            let value = if extended && flags & FLG_NONUT == 0 {
                crate::sidereal::true_ayanamsha_deg(mean, tt)?
            } else {
                mean
            };
            Ok((value, warning))
        };
        if !needs_anchor(selection) {
            return render(&SealedAyanamshaServing::new(selection, delta, None));
        }
        if let Some(rows) = rows {
            if self.ctx.mode() == crate::context::ProviderMode::SealedLeb {
                return crate::positions::star_serving::with_context_live_stored_anchor_provider(
                    &self.ctx,
                    self.delta_t_model(),
                    self.ctx.effective_iers_delta_t(),
                    rows,
                    |anchors| {
                        render(&SealedAyanamshaServing::new(
                            selection,
                            delta,
                            Some(anchors),
                        ))
                    },
                );
            }
        }
        // The defining computation helper selects the vector LEB adapter only
        // in forced LEB mode. An Auto live cache miss selects the main native
        // ephemeris, independently of a reader retained for its mean-UT prefix.
        // This cheap configuration clone admits no image; a cache hit remains
        // valid after reset/close without acquiring main. Standalone callers
        // with no live supplier retain their existing provider selection.
        let anchor_context =
            (rows.is_some() && self.ctx.mode() == crate::context::ProviderMode::Auto).then(|| {
                self.ctx
                    .clone()
                    .with_mode(crate::context::ProviderMode::NativeJpl)
            });
        let context = anchor_context.as_ref().unwrap_or(&self.ctx);
        crate::positions::star_serving::with_context_star_serving_captured_with_main(
            context,
            &self.star_main,
            None,
            self.delta_t_model(),
            self.ctx.effective_iers_delta_t(),
            0,
            None,
            |serving, _| {
                let anchors = serving.anchor_provider();
                let anchors = match rows {
                    Some(rows) => anchors.with_live_rows(rows),
                    None => anchors,
                };
                render(&SealedAyanamshaServing::new(
                    selection,
                    delta,
                    Some(&anchors),
                ))
            },
        )
        .result
    }

    /// Attempt the plain mean UT1 entry using the actual bound reader view.
    /// The stored table is queried at original UT1,
    /// then added once. Dynamic/ECL_DATE helpers invert that TT with configured
    /// Delta-T before reacquiring independent observation TT. USER_UT is skipped
    /// by the caller before this attempt; the local USER anchor remains TT.
    /// Reader assembly and resource errors retain their original ownership.
    fn public_fast_mean_ut(
        &self,
        ut: f64,
        rows: Option<&dyn AnchorRowAcquisition>,
    ) -> Result<Option<(f64, Option<Warning>)>> {
        if matches!(
            self.ctx.mode(),
            crate::context::ProviderMode::NativeJpl | crate::context::ProviderMode::HorizonsRemote
        ) || !self.ctx.leb_files().iter().any(|file| !file.is_closed())
        {
            return Ok(None);
        }
        let (selection, _) = crate::houses::sidereal_selection_for_context(&self.ctx);
        // Reader ownership is independent of the configured Native/Horizons
        // computation mode. This view changes only stack admission, never the
        // retained session or the delegated anchor provider's mode.
        let reader_context = self
            .ctx
            .clone()
            .with_mode(crate::context::ProviderMode::SealedLeb);
        crate::positions::context_serving::with_sealed_serving_stack(
            &reader_context,
            self.delta_t_model(),
            self.ctx.effective_iers_delta_t(),
            true,
            |stack| {
                let attempt = (|| {
                    let tt = ut + stack.serving.provider().delta_t(ut)?;
                    let delegated = DYNAMIC_AYANAMSHA_MODES.contains(&selection.mode())
                        || selection.bits() & SIDBIT_ECL_DATE != 0;
                    if delegated {
                        let helper_ut = tt - self.deltat(tt)?;
                        let observation_tt = crate::time::default_timescale_tt_from_ut1(helper_ut)?;
                        self.public_ayanamsha_at_tt(observation_tt, 0, false, rows)
                    } else {
                        self.public_ayanamsha_at_tt(tt, 0, false, rows)
                    }
                })();
                match attempt {
                    Ok(value) => Ok(Some(value)),
                    Err(error) => {
                        // A retained catalogue owner can classify its own acquisition error;
                        // otherwise use the native reader-miss categories.
                        let caught = match rows {
                            Some(rows) => rows.fast_mean_miss(&error)?,
                            None => None,
                        }
                        .unwrap_or_else(|| fast_miss(&error));
                        if caught { Ok(None) } else { Err(error) }
                    }
                }
            },
        )
    }

    /// Extended mean/true offset after the caller acquired the helper UT1.
    ///
    /// `helper_ut` is the actual independent TT→UT1 result of the released
    /// `_calc_ayanamsa_ex_value` prefix. It is not the original public UT1
    /// argument: an extended UT caller first constructs TT independently.
    /// This entry consumes only the remaining independent observation TT and
    /// selected geometry. It does not invert TT again. `flags` supplies the
    /// consumed no-aberration/no-deflection/NONUT bits; bindings keep the raw
    /// caller flags for their later exact return-flag operations. Session/C
    /// callers can pass their ordinary full flag word. Image/clock/error scope
    /// follows [`Self::public_ayanamsha`].
    pub fn public_ayanamsha_from_helper_ut(&self, helper_ut: f64, flags: i32) -> Result<Ayanamsha> {
        self.public_ayanamsha_from_helper_ut_with_anchor_rows(helper_ut, flags, None)
    }

    /// Evaluate an already acquired helper UT1 with optional live anchor rows.
    ///
    /// Epochs are Julian days, angles are degrees and flag semantics follow
    /// [`Self::public_ayanamsha_from_helper_ut`]. The callback is borrowed only
    /// through this synchronous calculation and is consumed solely by a live
    /// stellar/galactic mode, after independent observation TT acquisition.
    /// Classical, USER and formula modes do not read it. Callback failures
    /// propagate without a baked catalog substitute; its cache and operand
    /// ownership remain with the caller. `None` preserves static SDK behavior.
    pub fn public_ayanamsha_from_helper_ut_with_anchor_rows(
        &self,
        helper_ut: f64,
        flags: i32,
        rows: Option<&dyn AnchorRowAcquisition>,
    ) -> Result<Ayanamsha> {
        let (selection, _) = crate::houses::sidereal_selection_for_context(&self.ctx);
        self.public_ayanamsha_selected_from_helper_ut(helper_ut, flags, selection, rows)
    }

    /// Extended offset with an outer mode retained before mutable flag calls.
    ///
    /// The released module wrapper acquires `outer_mode` before interpreting
    /// raw flags. Its mean helper subsequently reads the current projection
    /// bits and USER reference epoch/value. This session must therefore be
    /// acquired after those flag calls: it supplies those late state words,
    /// configured defining-epoch clocks and actual anchor image ownership.
    /// `helper_ut` is already the independently acquired helper UT1, so this
    /// entry performs only the remaining UT1→TT observation conversion. The
    /// ordinary frozen Session/C entry retains its original full selection.
    pub fn public_ayanamsha_from_helper_ut_with_mode(
        &self,
        helper_ut: f64,
        flags: i32,
        outer_mode: i32,
    ) -> Result<Ayanamsha> {
        self.public_ayanamsha_from_helper_ut_with_mode_with_anchor_rows(
            helper_ut, flags, outer_mode, None,
        )
    }

    /// Retain an outer mode while consuming the current live anchor catalog.
    ///
    /// Current projection bits and USER epoch/value come from this session;
    /// `outer_mode` remains the earlier caller selection. Helper UT1, flags,
    /// output units and errors follow
    /// [`Self::public_ayanamsha_from_helper_ut_with_anchor_rows`]. Holding the
    /// callback introduces no field, cache, provider or clock read before the
    /// selected dynamic geometry. Its borrow does not escape this call.
    pub fn public_ayanamsha_from_helper_ut_with_mode_with_anchor_rows(
        &self,
        helper_ut: f64,
        flags: i32,
        outer_mode: i32,
        rows: Option<&dyn AnchorRowAcquisition>,
    ) -> Result<Ayanamsha> {
        let (current, _) = crate::houses::sidereal_selection_for_context(&self.ctx);
        let (_, t0, ayan_t0) = current.full();
        let mut selection = SiderealSelection::unset();
        selection.set_sidereal_mode(outer_mode, t0, ayan_t0);
        selection.set_retained_bits(current.bits());
        self.public_ayanamsha_selected_from_helper_ut(helper_ut, flags, selection, rows)
    }

    /// Complete the selected mean/true geometry after helper UT1 acquisition.
    /// Keeping selection explicit separates mutable module ownership from a
    /// frozen SDK session without altering either observation clock.
    fn public_ayanamsha_selected_from_helper_ut(
        &self,
        helper_ut: f64,
        flags: i32,
        selection: SiderealSelection,
        rows: Option<&dyn AnchorRowAcquisition>,
    ) -> Result<Ayanamsha> {
        let observation_tt = crate::time::default_timescale_tt_from_ut1(helper_ut)?;
        let (value, warning) =
            self.public_ayanamsha_selected_at_tt(observation_tt, flags, true, selection, rows)?;
        Ok(Ayanamsha {
            value,
            returned_flags: ayanamsha_ex_retflag(flags, selection.mode()),
            diagnostic: warning.map_or_else(String::new, |warning| warning.to_string()),
        })
    }

    /// Evaluate mean or true ayanamsha on the selected public time scale.
    ///
    /// `is_ut` selects UT1 or TT input; `extended` selects the returned flags and
    /// true-ayanamsha nutation term. Clock conversions follow the configured
    /// entry route. Dynamic stellar/galactic anchors use retained catalogue
    /// resources; unavailable anchors produce typed errors. Reader misses that
    /// permit analytical evaluation do not substitute fixed anchor offsets.
    pub fn public_ayanamsha(
        &self,
        epoch: f64,
        flags: i32,
        is_ut: bool,
        extended: bool,
    ) -> Result<Ayanamsha> {
        self.public_ayanamsha_with_anchor_rows(epoch, flags, is_ut, extended, None)
    }

    /// Plain or extended scalar ayanamsha with caller-owned live catalog rows.
    ///
    /// `epoch` is TT or UT1 Julian day according to `is_ut`; `extended`
    /// selects mean/true and returned-flag semantics. TT inputs reach the
    /// geometry directly; UT1 entries acquire their observation TT first, as
    /// in [`Self::public_ayanamsha`]. No row is read for catalog-free modes
    /// or unused defining epochs. A supplied row or cache failure remains the
    /// actual result; it does not authorize static catalog fallback. Immutable
    /// session resources survive the entire synchronous borrowed calculation.
    pub fn public_ayanamsha_with_anchor_rows(
        &self,
        epoch: f64,
        flags: i32,
        is_ut: bool,
        extended: bool,
        rows: Option<&dyn AnchorRowAcquisition>,
    ) -> Result<Ayanamsha> {
        let (selection, _) = crate::houses::sidereal_selection_for_context(&self.ctx);
        let (value, warning) = if !is_ut {
            // TT is already the observation time. Converting it through UT1
            // would let a reader's Delta-T table or another clock move the
            // supplied epoch, including an epoch-pair's defining instant.
            self.public_ayanamsha_at_tt(epoch, flags, extended, rows)?
        } else if extended {
            let entry_tt = crate::time::default_timescale_tt_from_ut1(epoch)?;
            let helper_ut = crate::time::default_timescale_ut1_from_tt(entry_tt)?;
            return self.public_ayanamsha_from_helper_ut_with_anchor_rows(helper_ut, flags, rows);
        } else {
            let ut = epoch;
            let user_ut = selection.mode() == SIDM_USER && selection.bits() & SIDBIT_USER_UT != 0;
            match if user_ut {
                None
            } else {
                self.public_fast_mean_ut(ut, rows)?
            } {
                Some(value) => value,
                None => {
                    let tt = crate::time::default_timescale_tt_from_ut1(ut)?;
                    self.public_ayanamsha_at_tt(tt, 0, false, rows)?
                }
            }
        };
        Ok(Ayanamsha {
            value,
            returned_flags: if extended {
                ayanamsha_ex_retflag(flags, selection.mode())
            } else {
                0
            },
            diagnostic: warning.map_or_else(String::new, |warning| warning.to_string()),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::RefCell;

    #[test]
    fn explicit_tt_preserves_every_defining_pair_with_or_without_reader_clocks() {
        let empty = crate::context::EphemerisContext::new();
        let stored = historical_reader_session().context().clone();
        for mut context in [empty, stored] {
            // An unrelated configured UT clock must never move an explicit TT
            // epoch. The reader also contains its own independent clock table.
            context.set_user_delta_t(Some(4.0));
            for (mode, offset, tt) in crate::sidereal::AYANAMSHA_DEFINING {
                context.set_sidereal_mode(mode, 0., 0.);
                let session = Session::from_context(context.clone());
                for extended in [false, true] {
                    let value = session
                        .public_ayanamsha(tt, FLG_NONUT, false, extended)
                        .unwrap();
                    let difference = (value.value - offset + 180.).rem_euclid(360.) - 180.;
                    assert!(
                        difference.abs() * 3600. <= 0.001,
                        "mode={mode}, TT={tt}: {difference}"
                    );
                }
            }
        }
    }

    #[test]
    fn acquired_mean_tt_consumes_only_selected_defining_epoch_clock() {
        let session = Session::from_context(crate::context::EphemerisContext::new());
        let mut selection = SiderealSelection::unset();
        selection.set_sidereal_mode(0, 2451545., 0.);
        assert!(
            session
                .public_mean_at_observation_tt_with_selection(
                    2451545.,
                    selection,
                    None,
                    &|_| panic!("classical defining epoch has no model clock"),
                )
                .is_ok()
        );
        selection.set_sidereal_mode(SIDM_VALENS_MOON, 2451545., 0.);
        let seen = RefCell::new(Vec::new());
        let marker = Error::MathDomain { value: 12.5 };
        let result = session.public_mean_at_observation_tt_with_selection(
            2451545.,
            selection,
            None,
            &|ut| {
                seen.borrow_mut().push(ut);
                Err(marker.clone())
            },
        );
        assert_eq!(result.unwrap_err(), marker);
        assert_eq!(*seen.borrow(), [crate::sidereal::VALENS_MOON_T0_UT]);
    }
    use crate::stars::{AcquiredAnchorRow, AnchorSample, AnchorSampleOptions, NativeStarInputs};
    use std::sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    };

    struct CachedRow;
    impl AcquiredAnchorRow for CachedRow {
        fn sample(&self, _tt: f64, _options: AnchorSampleOptions) -> Result<AnchorSample> {
            Ok(AnchorSample::Cached(195.0))
        }
        fn galactic_inputs(&self) -> Result<NativeStarInputs> {
            panic!("Spica longitude must not consume Galactic Center keywords")
        }
    }

    struct RetryRows {
        lookups: AtomicUsize,
        catches: AtomicUsize,
        catch: bool,
    }

    struct CacheOnlyRows(AtomicUsize);
    impl AnchorRowAcquisition for CacheOnlyRows {
        fn acquire_row(&self, key: &str) -> Result<Arc<dyn AcquiredAnchorRow>> {
            assert_eq!(key, "SPICA");
            self.0.fetch_add(1, Ordering::SeqCst);
            Ok(Arc::new(CachedRow))
        }
    }

    #[test]
    fn live_cached_extended_anchor_survives_without_a_stored_main_image() {
        for mode in [
            crate::context::ProviderMode::SealedLeb,
            crate::context::ProviderMode::Auto,
        ] {
            let mut context = crate::context::EphemerisContext::new()
                .with_mode(mode)
                .with_leb_file(None);
            assert!(context.leb_files().is_empty());
            context.set_sidereal_mode(crate::constants::SIDM_TRUE_CITRA, 0.0, 0.0);
            let session = Session::from_context(context);
            let rows = CacheOnlyRows(AtomicUsize::new(0));
            for is_ut in [false, true] {
                let value = session
                    .public_ayanamsha_with_anchor_rows(
                        2460400.25,
                        FLG_NONUT,
                        is_ut,
                        true,
                        Some(&rows),
                    )
                    .expect("complete cached anchor must not acquire an absent main image");
                assert_eq!(value.value, 15.0);
            }
            assert_eq!(rows.0.load(Ordering::SeqCst), 2);
        }
    }
    impl AnchorRowAcquisition for RetryRows {
        fn acquire_row(&self, key: &str) -> Result<Arc<dyn AcquiredAnchorRow>> {
            assert_eq!(key, "SPICA");
            if self.lookups.fetch_add(1, Ordering::SeqCst) == 0 {
                Err(Error::ConfigurationMissing {
                    missing_config: "synthetic_live_anchor_row_failure",
                    suggestion: "classify the retained same-call row refusal",
                })
            } else {
                Ok(Arc::new(CachedRow))
            }
        }
        fn fast_mean_miss(&self, error: &Error) -> Result<Option<bool>> {
            assert!(matches!(
                error.primary_error(),
                Error::ConfigurationMissing {
                    missing_config: "synthetic_live_anchor_row_failure",
                    ..
                }
            ));
            self.catches.fetch_add(1, Ordering::SeqCst);
            Ok(Some(self.catch))
        }
    }

    fn historical_reader_session() -> Session {
        use sha2::{Digest, Sha256};
        let path = std::env::var_os("LIBEPHEMERIS_TEST_LEB2")
            .map(std::path::PathBuf::from)
            .unwrap_or_else(crate::test_data::core_leb_path);
        let raw = std::fs::read(&path).expect("historical fixture is required, not skipped");
        let digest = Sha256::digest(&raw)
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect::<String>();
        assert!(crate::test_data::admitted_core_hash(&(digest)));
        let mut context = crate::context::EphemerisContext::new()
            .with_mode(crate::context::ProviderMode::SealedLeb);
        context.attach_leb_file(crate::leb::LebFile::open(&path).unwrap());
        context.set_sidereal_mode(crate::constants::SIDM_TRUE_CITRA, 0.0, 0.0);
        Session::from_context(context)
    }

    #[test]
    fn live_fast_mean_catch_retries_the_same_supplier_without_static_substitution() {
        let session = historical_reader_session();
        let rows = RetryRows {
            lookups: AtomicUsize::new(0),
            catches: AtomicUsize::new(0),
            catch: true,
        };
        let value = session
            .public_ayanamsha_with_anchor_rows(2460400.25, 0, true, false, Some(&rows))
            .unwrap();
        // Constructed cache longitude195 minus the Citra target180 is15;
        // this is an algebraic controller test, not native-derived science.
        assert_eq!(value.value, 15.0);
        assert_eq!(rows.lookups.load(Ordering::SeqCst), 2);
        assert_eq!(rows.catches.load(Ordering::SeqCst), 1);
        assert_eq!(
            session.context().sidereal_mode(),
            crate::constants::SIDM_TRUE_CITRA
        );
    }

    #[test]
    fn live_fast_mean_unhandled_refusal_remains_terminal() {
        let session = historical_reader_session();
        let rows = RetryRows {
            lookups: AtomicUsize::new(0),
            catches: AtomicUsize::new(0),
            catch: false,
        };
        let error = session
            .public_ayanamsha_with_anchor_rows(2460400.25, 0, true, false, Some(&rows))
            .unwrap_err();
        assert!(matches!(
            error.primary_error(),
            Error::ConfigurationMissing {
                missing_config: "synthetic_live_anchor_row_failure",
                ..
            }
        ));
        assert_eq!(rows.lookups.load(Ordering::SeqCst), 1);
        assert_eq!(rows.catches.load(Ordering::SeqCst), 1);
    }

    #[test]
    fn live_catalog_owner_is_unused_by_classical_and_formula_modes() {
        struct Unused;
        impl AnchorRowAcquisition for Unused {
            fn acquire_row(&self, _key: &str) -> Result<Arc<dyn AcquiredAnchorRow>> {
                panic!("catalog-free mode must not read a row")
            }
        }
        for mode in [0, SIDM_GALALIGN_MARDYKS, SIDM_VALENS_MOON, SIDM_USER, -1] {
            let mut context = crate::context::EphemerisContext::new();
            context.set_sidereal_mode(mode, 2451545.0, 9.0);
            let session = Session::from_context(context);
            for (is_ut, extended) in [(false, false), (true, false), (false, true), (true, true)] {
                let live = session
                    .public_ayanamsha_with_anchor_rows(
                        2460400.25,
                        FLG_NONUT,
                        is_ut,
                        extended,
                        Some(&Unused),
                    )
                    .unwrap();
                let ordinary = session
                    .public_ayanamsha(2460400.25, FLG_NONUT, is_ut, extended)
                    .unwrap();
                assert_eq!(live.value.to_bits(), ordinary.value.to_bits());
                assert_eq!(live.returned_flags, ordinary.returned_flags);
                assert_eq!(live.diagnostic, ordinary.diagnostic);
            }
        }
    }

    #[test]
    fn empty_auto_classical_offsets_need_no_reader_or_main() {
        let session = Session::from_context(crate::context::EphemerisContext::new());
        for (is_ut, extended) in [(true, false), (false, false), (true, true), (false, true)] {
            let result = session
                .public_ayanamsha(2451545.0, FLG_NONUT, is_ut, extended)
                .unwrap();
            assert!(result.value.is_finite());
            assert!(session.context().spk_registry().is_empty());
            assert!(session.context().leb_files().is_empty());
        }
    }

    #[test]
    fn model_override_does_not_replace_independent_observation_clock() {
        let ordinary = Session::from_context(crate::context::EphemerisContext::new());
        let mut overridden = ordinary.context().clone();
        overridden.set_user_delta_t(Some(4.0));
        let overridden = Session::from_context(overridden);
        for (is_ut, extended) in [(true, false), (false, false), (true, true), (false, true)] {
            let a = ordinary
                .public_ayanamsha(2460400.25, 0, is_ut, extended)
                .unwrap();
            let b = overridden
                .public_ayanamsha(2460400.25, 0, is_ut, extended)
                .unwrap();
            assert_eq!(a.value.to_bits(), b.value.to_bits());
        }
    }

    #[test]
    fn mutable_outer_user_mode_consumes_late_pair_and_bits() {
        let helper_ut = 2460000.0;
        let tt = crate::time::default_timescale_tt_from_ut1(helper_ut).unwrap();
        let mut context = crate::context::EphemerisContext::new();
        // A caller changed the current base mode as well as the USER pair.
        // The outer public call still selects USER; at its defining TT the
        // Method-B increment and ECL_DATE projection are exactly zero.
        context.set_sidereal_mode(crate::constants::SIDM_J2000 | SIDBIT_ECL_DATE, tt, 9.0);
        let session = Session::from_context(context);
        let mutable = session
            .public_ayanamsha_from_helper_ut_with_mode(helper_ut, FLG_NONUT, SIDM_USER)
            .unwrap();
        assert_eq!(mutable.value, 9.0);
        assert_eq!(mutable.returned_flags, 66);
        let frozen = session
            .public_ayanamsha_from_helper_ut(helper_ut, FLG_NONUT)
            .unwrap();
        assert_ne!(frozen.value.to_bits(), mutable.value.to_bits());
        assert_eq!(
            session.context().sidereal_mode(),
            crate::constants::SIDM_J2000
        );
    }

    #[test]
    fn mutable_outer_user_mode_consumes_late_user_ut_model() {
        let helper_ut = 2460000.0;
        let tt = crate::time::default_timescale_tt_from_ut1(helper_ut).unwrap();
        let mut context = crate::context::EphemerisContext::new();
        context.set_user_delta_t(Some(4.0));
        context.set_sidereal_mode(crate::constants::SIDM_J2000 | SIDBIT_USER_UT, tt - 4.0, 9.0);
        let session = Session::from_context(context);
        let result = session
            .public_ayanamsha_from_helper_ut_with_mode(helper_ut, FLG_NONUT, SIDM_USER)
            .unwrap();
        assert_eq!(result.value, 9.0);
        let mut without_bits = session.context().clone();
        without_bits.set_sidereal_bits(0);
        let without_bits = Session::from_context(without_bits)
            .public_ayanamsha_from_helper_ut_with_mode(helper_ut, FLG_NONUT, SIDM_USER)
            .unwrap();
        assert_ne!(result.value.to_bits(), without_bits.value.to_bits());
    }
}
