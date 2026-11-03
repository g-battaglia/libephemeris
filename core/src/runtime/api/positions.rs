//! Full position entry composition; contexts, clocks and providers stay owned.
use super::*;
use crate::houses::sidereal_selection_for_context;
use crate::positions::context_serving::{self as serving, ServingStack};
use crate::positions::{self, TopoPlan};
use crate::sidereal::{
    SIDEREAL_SPEED_HALF_STEP_DAYS, SealedAyanamshaServing, SiderealRoute,
    apply_sidereal_correction, resolve_sidereal_route, strip_to_base_mode,
};

#[path = "position_acquisitions.rs"]
mod acquisition;
pub use acquisition::PositionAcquisitions;
use acquisition::{FastPositionAnchors, PositionAnchors};

#[path = "position_continuation.rs"]
mod continuation;
#[path = "pctr_main_image.rs"]
mod pctr_main_image;
#[path = "position_request.rs"]
mod request;
use request::{InputEpoch, PositionRequest};

#[path = "position_catalog.rs"]
mod catalog;
#[path = "position_dispatch.rs"]
mod dispatch;
#[path = "position_stored_attempt.rs"]
mod stored_attempt;

/// Unshaped spherical model callback consumed by the public position entry.
type ModelCompute<'a> = dyn Fn(f64, i32, i32) -> Result<([f64; 6], i32), Error> + 'a;

/// Lazy UT1 partner for one stored-state TT site sample, in Julian days.
type SiteClock<'a> = dyn Fn(f64) -> Result<f64, Error> + 'a;

#[cfg(test)]
std::thread_local! {
    /// Scoped observation of dispatch only; it never replaces a source.
    static OWNED_MINOR_ENTRY_PROBE: std::cell::RefCell<Option<Vec<bool>>> = const {
        std::cell::RefCell::new(None)
    };
}

impl Session {
    #[cfg(test)]
    pub(crate) fn observe_minor_entry_ownership<T>(
        operation: impl FnOnce() -> T,
    ) -> (T, Vec<bool>) {
        struct Reset;
        impl Drop for Reset {
            fn drop(&mut self) {
                OWNED_MINOR_ENTRY_PROBE.with(|probe| *probe.borrow_mut() = None);
            }
        }
        OWNED_MINOR_ENTRY_PROBE.with(|probe| {
            assert!(probe.borrow().is_none(), "minor-entry probe already active");
            *probe.borrow_mut() = Some(Vec::new());
        });
        let _reset = Reset;
        let result = operation();
        let observations = OWNED_MINOR_ENTRY_PROBE.with(|probe| probe.borrow_mut().take().unwrap());
        (result, observations)
    }

    /// Calculate a position from a UT1 Julian day.
    ///
    /// # Arguments
    ///
    /// - `jd_ut`: UT1 Julian day, not UTC or TT. Ordinary stored-body requests
    ///   convert it using this session's model Delta-T. Native observer/frame,
    ///   stellar and `ECL_NUT` calculations retain their independently defined
    ///   clocks; a user Delta-T override does not replace every such clock.
    /// - `body`: a public identifier from [`crate::constants`]. Families include
    ///   planets, lunar points, minor/hypothetical bodies, moons, catalog stars,
    ///   angles and cached Arabic parts. Negative mean/true-node identifiers
    ///   request the descending node; recognized `AST_OFFSET + N` aliases are
    ///   remapped. Actual support and coverage depend on the selected data/model.
    /// - `flags`: bitwise OR of `FLG_*` constants. Zero requests an apparent
    ///   geocentric ecliptic place of date without rates. `FLG_EQUATORIAL`
    ///   selects equatorial coordinates, `FLG_XYZ` Cartesian coordinates and
    ///   `FLG_RADIANS` spherical radians. J2000, no-nutation and sidereal flags
    ///   select their corresponding frames/projections. Topocentric output
    ///   requires a configured observer when that body consumes a surface site.
    ///
    /// For ordinary positions, center priority is TOPOCTR > BARYCTR > HELCTR >
    /// geocentric. TRUEPOS omits light time, aberration and deflection; NOABERR
    /// and NOGDEFL suppress their respective corrections. Ephemeris label bits
    /// do not override the context's provider mode. `FLG_SPEED3` alone requests
    /// rates; explicit `FLG_SPEED` supersedes it.
    ///
    /// # Returns
    ///
    /// An owned [`Position`], whose documentation specifies component indices,
    /// units, returned flags and the `ECL_NUT` exception. Position rates remain
    /// **per TT day**, not per UT1 day. Angle/Arabic-part identifiers return a
    /// longitude in degrees followed by five zeros rather than ordinary body
    /// output shaping; Arabic parts consume the context's angle cache, not the
    /// requested epoch. Catalog-star calculations follow [`Self::fixstar_ut`].
    ///
    /// # Errors
    ///
    /// Returns typed errors for nonfinite consumed inputs, unsupported flag
    /// combinations (including an exclusively MOSEPH barycentric request),
    /// unsupported bodies, missing consumed observer/cache/data, actual
    /// provider/model coverage misses, and clock/frame/numerical failures.
    /// File or permitted network acquisition can fail when the selected route
    /// needs resources. Sealed LEB mode does not silently replace a missing
    /// stored planetary channel with JPL. Checks occur at the consuming branch;
    /// cached parts and other data-independent results need not consume time.
    ///
    /// # Example
    ///
    /// This explicitly opens an existing local file; it never downloads one.
    ///
    /// ```no_run
    /// use libephemeris_core::{
    ///     api,
    ///     context::{EphemerisContext, ProviderMode},
    ///     leb::LebFile,
    /// };
    ///
    /// fn main() -> Result<(), Box<dyn std::error::Error>> {
    ///     let mut context = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
    ///     let file = LebFile::open(std::path::Path::new("/path/to/base_core.leb2"))?;
    ///     context.attach_leb_file(file);
    ///     let session = context.primary_session_builder().build()?;
    ///     let sun = session.calc_ut(2451545.0, api::SUN, api::FLG_SPEED)?;
    ///     println!("Longitude: {} degrees", sun.longitude());
    ///     println!("Longitude rate: {} degrees/TT day", sun.values[3]);
    ///     Ok(())
    /// }
    /// ```
    pub fn calc_ut(&self, jd_ut: f64, body: i32, flags: i32) -> Result<Position, Error> {
        self.position_entry(InputEpoch::Ut1(jd_ut), body, flags)
    }
    /// Calculate a position from a TT Julian day.
    ///
    /// `jd_tt` is Terrestrial Time in Julian days; no initial UT1-to-TT model
    /// conversion is applied. Body identifiers, correction/frame/representation
    /// flags, resource requirements and errors follow [`Self::calc_ut`].
    /// Independent observer/frame clocks and sidereal defining epochs retain
    /// their own conversions. For `ECL_NUT`, the supplied TT directly selects
    /// the frame angles, with no body-position shaping.
    ///
    /// Returns an owned [`Position`] with ordinary rates per TT day. Unlike
    /// `calc_ut`, the TT entry returns no ephemeris-selection bit when none was
    /// requested; returned labels still do not identify the physical provider.
    pub fn calc(&self, jd_tt: f64, body: i32, flags: i32) -> Result<Position, Error> {
        self.position_entry(InputEpoch::Tt(jd_tt), body, flags)
    }
    pub(crate) fn calc_from_native_time_scalar(
        &self,
        jd: f64,
        body: i32,
        flags: i32,
    ) -> Result<Position, Error> {
        self.calc(jd, body, flags)
    }
    fn position_entry(
        &self,
        epoch: InputEpoch,
        body: i32,
        input_flags: i32,
    ) -> Result<Position, Error> {
        let mut supplier = super::NativeMinorAssist::discover();
        let acquisitions = PositionAcquisitions {
            minor_assist: Some(&supplier),
            ..Default::default()
        };
        let mut position = self.position_entry_with_owned_minor(
            &acquisitions,
            None,
            Some(&supplier),
            PositionRequest {
                epoch,
                body,
                input_flags,
            },
        )?;
        let diagnostics = supplier.drain_diagnostics_exclusive();
        position.warnings.extend(
            diagnostics
                .warnings
                .into_iter()
                .map(super::Advisory::from_engine),
        );
        position.messages.extend(diagnostics.messages);
        Ok(position)
    }
    /// Calculate a TT position with optional, synchronously borrowed providers.
    ///
    /// Epoch, body, flags, layout, units and errors follow [`Self::calc`].
    /// `acquisitions` supplies only the selected model, live sidereal catalog
    /// or numerical minor-body stages described by [`PositionAcquisitions`].
    /// Unused providers are not called. The caller retains them until return;
    /// consumed provider errors propagate without assuming callback purity.
    pub fn calc_with_acquisitions(
        &self,
        jd: f64,
        body: i32,
        flags: i32,
        acquisitions: &PositionAcquisitions<'_>,
    ) -> Result<Position, Error> {
        self.position_entry_with_owned_minor(
            acquisitions,
            None,
            None,
            PositionRequest {
                epoch: InputEpoch::Tt(jd),
                body,
                input_flags: flags,
            },
        )
    }
    /// UT1 counterpart of [`Self::calc_with_acquisitions`].
    ///
    /// `jd` is a UT1 Julian day. Position units, rates per TT day, flag handling
    /// and errors follow [`Self::calc_ut`]; provider borrowing and laziness
    /// follow [`Self::calc_with_acquisitions`].
    pub fn calc_ut_with_acquisitions(
        &self,
        jd: f64,
        body: i32,
        flags: i32,
        acquisitions: &PositionAcquisitions<'_>,
    ) -> Result<Position, Error> {
        self.position_entry_with_owned_minor(
            acquisitions,
            None,
            None,
            PositionRequest {
                epoch: InputEpoch::Ut1(jd),
                body,
                input_flags: flags,
            },
        )
    }
    /// Calculate a TT position with a synchronously borrowed model callback.
    ///
    /// Arguments, return units/layout and errors follow [`Self::calc`].
    /// `models` follows [`positions::AnalyticalModelCallback`]; its numeric
    /// result is owned, but the callback itself is borrowed only until return.
    /// It is called only by a selected analytical-model stage, possibly at
    /// multiple rate/retarded samples. `Ok(None)` keeps that stage's built-in
    /// model, while `Err` terminates the consumed sample.
    pub fn calc_with_models(
        &self,
        jd: f64,
        body: i32,
        flags: i32,
        models: &positions::AnalyticalModelCallback<'_>,
    ) -> Result<Position, Error> {
        self.calc_with_acquisitions(
            jd,
            body,
            flags,
            &PositionAcquisitions {
                models: Some(models),
                ..Default::default()
            },
        )
    }
    /// UT1 counterpart of [`Self::calc_with_models`].
    ///
    /// `jd` is a UT1 Julian day; output and errors follow [`Self::calc_ut`].
    /// Model callback ownership, repeated samples and optional-state behavior
    /// follow [`Self::calc_with_models`].
    pub fn calc_ut_with_models(
        &self,
        jd: f64,
        body: i32,
        flags: i32,
        models: &positions::AnalyticalModelCallback<'_>,
    ) -> Result<Position, Error> {
        self.calc_ut_with_acquisitions(
            jd,
            body,
            flags,
            &PositionAcquisitions {
                models: Some(models),
                ..Default::default()
            },
        )
    }
    /// Complete UT/TT dispatch; one temporary serving view over shared data.
    pub(super) fn position_entry_with_acquisitions_and_main_image(
        &self,
        acquisitions: &PositionAcquisitions<'_>,
        main_image: Option<&dyn positions::MainImageAcquisition>,
        jd: f64,
        body: i32,
        flags: i32,
        is_ut: bool,
    ) -> Result<Position, Error> {
        self.position_entry_with_owned_minor(
            acquisitions,
            main_image,
            None,
            PositionRequest {
                epoch: InputEpoch::from_entry(jd, is_ut),
                body,
                input_flags: flags,
            },
        )
    }

    /// Runtime's concrete numerical owner survives the late main-image join.
    /// Public acquisition entry points continue to treat all suppliers as opaque.
    pub(crate) fn position_with_native_minor_and_main_image(
        &self,
        supplier: &super::NativeMinorAssist<'_>,
        main_image: &dyn positions::MainImageAcquisition,
        jd: f64,
        body: i32,
        flags: i32,
        is_ut: bool,
    ) -> Result<Position, Error> {
        let acquisitions = PositionAcquisitions {
            minor_assist: Some(supplier),
            ..Default::default()
        };
        self.position_entry_with_owned_minor(
            &acquisitions,
            Some(main_image),
            Some(supplier),
            PositionRequest {
                epoch: InputEpoch::from_entry(jd, is_ut),
                body,
                input_flags: flags,
            },
        )
    }

    /// Concrete native ownership is a private default-entry fact. Public
    /// acquisition records remain opaque and retain their original call shape.
    fn position_entry_with_owned_minor(
        &self,
        acquisitions: &PositionAcquisitions<'_>,
        main_image: Option<&dyn positions::MainImageAcquisition>,
        owned_minor: Option<&super::NativeMinorAssist<'_>>,
        request: PositionRequest,
    ) -> Result<Position, Error> {
        #[cfg(test)]
        OWNED_MINOR_ENTRY_PROBE.with(|probe| {
            if let Some(observations) = probe.borrow_mut().as_mut() {
                observations.push(owned_minor.is_some());
            }
        });
        let PositionRequest {
            epoch,
            body,
            input_flags,
        } = request;
        let (jd, is_ut) = epoch.into_entry();
        // Catalog identifiers finish before any body or live-anchor serving.
        if let Some(position) = self.catalog_position(request)? {
            return Ok(position);
        }
        let owned_anchors = PositionAnchors {
            session: self,
            rows: acquisitions.anchor_rows,
        };
        let live = Some(&owned_anchors as &dyn crate::sidereal::DynamicAnchorProvider);
        // ECL_NUT owns independent frame/time callbacks and consumes no body
        // image, configured model, observer or mutable analytical model feed.
        if body == ECL_NUT {
            let native = self.ctx.mode() == ProviderMode::NativeJpl;
            let result = (|| {
                if native {
                    positions::validate_barycentric_moseph(
                        input_flags,
                        if is_ut { "calc_ut" } else { "calc" },
                    )?;
                }
                independent_nutation_entry(epoch, input_flags)
            })();
            let (values, returned_flags) = result?;
            return Ok(completed_position(
                values,
                returned_flags,
                self.diagnostics.join("; "),
            ));
        }
        if self.ctx.mode() != ProviderMode::NativeJpl && !jd.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd });
        }
        if matches!(
            self.ctx.mode(),
            ProviderMode::NativeJpl | ProviderMode::HorizonsRemote
        ) || (self.ctx.mode() == ProviderMode::Auto && self.ctx.leb_files().is_empty())
        {
            let (values, returned_flags) =
                self.position_after_fast(request, acquisitions, live, main_image, owned_minor)?;
            return Ok(completed_position(
                values,
                returned_flags,
                self.diagnostics.join("; "),
            ));
        }
        if let Some((values, returned_flags)) =
            self.specialized_position(request, acquisitions, live)?
        {
            return Ok(completed_position(
                values,
                returned_flags,
                self.diagnostics.join("; "),
            ));
        }
        // Ordinary fast requests select their observer before configured time.
        // Raw center/ephemeris conflicts keep the public flag-validation order.
        positions::validate_barycentric_moseph(
            input_flags,
            if is_ut { "calc_ut" } else { "calc" },
        )?;
        // This snapshot precedes the stored attempt and any Auto continuation.
        let diagnostic = self.diagnostics.join("; ");
        let attempt = self.ordinary_stored_position_attempt(
            request,
            acquisitions,
            live,
            main_image.is_some(),
        );
        let (values, returned_flags) = match attempt {
            Ok(value) => value,
            Err(error)
                if self.ctx.mode() == ProviderMode::Auto
                    && matches!(
                        error.primary_error(),
                        Error::LebBodyNotFound { .. }
                            | Error::LebOutOfRange { .. }
                            | Error::SealedChannelNotServed { .. }
                            | Error::EphemerisRange { .. }
                    ) =>
            {
                self.position_after_fast(request, acquisitions, live, main_image, owned_minor)?
            }
            Err(error)
                if self.ctx.mode() == ProviderMode::SealedLeb
                    && crate::planetary_moons::is_planetary_moon(body)
                    && matches!(
                        error.primary_error(),
                        Error::LebBodyNotFound { .. } | Error::SealedChannelNotServed { .. }
                    ) =>
            {
                return Err(Error::ConfigurationMissing {
                    missing_config: "jpl_source_access",
                    suggestion: "this known moon has no LEB channel; select a mode with satellite SPK access",
                });
            }
            Err(error) => return Err(error),
        };
        Ok(completed_position(values, returned_flags, diagnostic))
    }

    /// Planet-centric position at scalar TT over the context's explicit data.
    /// Auto/Horizons use the selected native main directly; SealedLeb keeps
    /// its stored vector path. No fast retry or Horizons-vector leg is used.
    /// Base sidereal, fixed-epoch and ECL_T0/SSY_PLANE frames use the stored
    /// context selection. Split SPEED/retardation clocks, observer precedence,
    /// units, output flags and errors follow [`positions::calc_pctr_in_context`].
    pub fn calc_pctr(
        &self,
        jd_tt: f64,
        body: i32,
        center: i32,
        flags: i32,
    ) -> Result<Position, Error> {
        self.calc_pctr_dispatch(jd_tt, body, center, flags, None, None)
    }
    /// Module TT PCTR with live catalog rows only at the consumed sidereal tail.
    /// Body and SPEED geometry precede the independent ayanamsha Time.
    pub fn calc_pctr_with_anchor_rows(
        &self,
        jd_tt: f64,
        body: i32,
        center: i32,
        flags: i32,
        rows: Option<&dyn crate::stars::AnchorRowAcquisition>,
    ) -> Result<Position, Error> {
        self.calc_pctr_dispatch(jd_tt, body, center, flags, rows, None)
    }
}

/// Package completed numerical output without changing values or flags.
///
/// The caller owns diagnostic evaluation timing: ordinary stored attempts take
/// their snapshot before serving; catalog, nutation, specialized and PCTR routes
/// take it after completion. Native minor advisories are drained separately by
/// the entry that owns that supplier, never by this packaging step.
fn completed_position(values: [f64; 6], returned_flags: i32, diagnostic: String) -> Position {
    Position {
        warnings: Vec::new(),
        messages: Vec::new(),
        values,
        returned_flags,
        diagnostic,
    }
}

/// Own a completed PCTR result and this Session's existing advisories.
/// It performs no clock, flag, provider or numerical work.
pub(super) fn pctr_position(session: &Session, values: [f64; 6], returned_flags: i32) -> Position {
    completed_position(values, returned_flags, session.diagnostics.join("; "))
}

#[cfg(test)]
mod completed_position_tests {
    use super::*;

    #[test]
    fn packaging_preserves_completed_words_flags_and_diagnostic_snapshot() {
        let values = [
            -0.0,
            0.0,
            f64::from_bits(0x7ff8_0000_0000_1234),
            f64::INFINITY,
            f64::NEG_INFINITY,
            -1.25,
        ];
        // These are completed flags, not a new request to normalize.
        let returned_flags = i32::MIN | FLG_SPEED3 | FLG_SIDEREAL;
        let mut diagnostics = vec![String::from("first"), String::from("second")];
        let snapshot = diagnostics.join("; ");
        diagnostics.push(String::from("later"));
        let position = completed_position(values, returned_flags, snapshot);
        assert_eq!(position.values.map(f64::to_bits), values.map(f64::to_bits));
        assert_eq!(position.returned_flags, returned_flags);
        assert_eq!(position.diagnostic, "first; second");
        assert!(position.warnings.is_empty());
        assert!(position.messages.is_empty());
    }

    #[test]
    fn pctr_packaging_uses_the_same_completed_result_contract() {
        let mut session = Session::from_context(EphemerisContext::new());
        session.diagnostics = vec![String::from("first"), String::from("second")];
        let values = [-0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
        let returned_flags = FLG_SPEED3 | FLG_SIDEREAL;
        let position = pctr_position(&session, values, returned_flags);
        assert_eq!(position.values.map(f64::to_bits), values.map(f64::to_bits));
        assert_eq!(position.returned_flags, returned_flags);
        assert_eq!(position.diagnostic, "first; second");
        assert!(position.warnings.is_empty());
        assert!(position.messages.is_empty());
    }
}

/// Data-independent ECL_NUT entry using the public flag/echo reduction.
/// The UT callback owns its bundled UT1-to-TT clock; TT uses its supplied word.
fn independent_nutation_entry(epoch: InputEpoch, flags: i32) -> Result<([f64; 6], i32), Error> {
    let unavailable_body = |_: f64, body, _: i32| Err(Error::UnknownBody { body_id: body });
    let engine = positions::CalcEngine {
        ut1_to_tt_jd: &crate::time::default_timescale_tt_from_ut1,
        compute_tt: &unavailable_body,
        nutation_obliquity_ut: &positions::ecl_nut_result_ut,
        nutation_obliquity_tt: &|tt, nut_flags| {
            positions::ecl_nut_result(tt, nut_flags, &crate::erfa::nut06a)
        },
        jpl_fallback_active: false,
    };
    match epoch {
        InputEpoch::Ut1(jd) => positions::calc_ut(&engine, jd, ECL_NUT, flags),
        InputEpoch::Tt(jd) => positions::calc(&engine, jd, ECL_NUT, flags),
    }
}

/// Trace the actual topocentric Earth bypass without labeling ordinary engines.
fn earth_origin_trace<R>(
    body: i32,
    flags: i32,
    operation: impl FnOnce() -> Result<R, Error>,
) -> Result<R, Error> {
    if body == EARTH && flags & FLG_TOPOCTR != 0 {
        positions::with_serving_trace_source(body, "Analytical", false, operation)
    } else {
        operation()
    }
}

#[cfg(test)]
mod independent_nutation_tests {
    use super::*;
    use crate::trace::TraceEvent;
    use std::cell::Cell;

    #[test]
    fn empty_session_nutation_uses_independent_tt_and_ut_clocks() {
        let jd = 2451545.125;
        let flags = FLG_JPLEPH | FLG_SPEED | FLG_SPEED3 | FLG_XYZ | FLG_RADIANS;
        for mode in [
            ProviderMode::SealedLeb,
            ProviderMode::Auto,
            ProviderMode::NativeJpl,
            ProviderMode::HorizonsRemote,
        ] {
            let mut ctx = EphemerisContext::new()
                .with_mode(mode)
                .with_delta_t_model(Some(crate::time::DeltaTModel::EspenakMeeus));
            // A model override would destroy the selected epoch if consumed.
            // ECL_NUT instead owns the bundled independent Time callback.
            ctx.set_user_delta_t(Some(100.0));
            let session = Session::from_context(ctx);
            for is_ut in [false, true] {
                let tt = if is_ut {
                    crate::time::default_timescale_tt_from_ut1(jd).unwrap()
                } else {
                    jd
                };
                let expected = positions::ecl_nut_result(
                    tt,
                    FLG_JPLEPH | FLG_SPEED | FLG_XYZ | FLG_RADIANS,
                    &crate::erfa::nut06a,
                )
                .unwrap();
                let (result, events) = crate::trace::capture_events(|| {
                    if is_ut {
                        session.calc_ut(jd, ECL_NUT, flags)
                    } else {
                        session.calc(jd, ECL_NUT, flags)
                    }
                });
                let result = result.unwrap();
                assert_eq!(
                    result.values.map(f64::to_bits),
                    expected.0.map(f64::to_bits)
                );
                assert_eq!(result.returned_flags, expected.1);
                assert_eq!(
                    events,
                    vec![TraceEvent::Record {
                        body: ECL_NUT,
                        source: "ERFA".into()
                    }]
                );
            }
        }
    }

    #[test]
    fn model_nutation_never_consumes_model_or_ephemeris_callbacks() {
        let session = SessionBuilder::isolated().build().unwrap();
        let calls = Cell::new(0);
        let models = |_, _, _| {
            calls.set(calls.get() + 1);
            Err(Error::UnknownBody { body_id: 913 })
        };
        for is_ut in [false, true] {
            let result = if is_ut {
                session.calc_ut_with_models(2451545.0, ECL_NUT, FLG_SWIEPH, &models)
            } else {
                session.calc_with_models(2451545.0, ECL_NUT, FLG_SWIEPH, &models)
            }
            .unwrap();
            assert_eq!(
                [result.values[4].to_bits(), result.values[5].to_bits()],
                [0, 0]
            );
            assert_eq!(result.returned_flags, FLG_SWIEPH);
        }
        assert_eq!(calls.get(), 0);
    }

    #[test]
    fn data_free_nutation_retains_existing_guard_order_and_failed_trace() {
        let flags = FLG_BARYCTR | FLG_MOSEPH;
        for mode in [ProviderMode::SealedLeb, ProviderMode::NativeJpl] {
            let session = Session::from_context(EphemerisContext::new().with_mode(mode));
            for is_ut in [false, true] {
                let entry = if is_ut { "calc_ut" } else { "calc" };
                let (result, events) = crate::trace::capture_events(|| {
                    if is_ut {
                        session.calc_ut(f64::NAN, ECL_NUT, flags)
                    } else {
                        session.calc(f64::NAN, ECL_NUT, flags)
                    }
                });
                if mode == ProviderMode::NativeJpl {
                    let error = result.unwrap_err();
                    assert_eq!(
                        error.primary_error(),
                        &Error::UnsupportedFlagCombination { entry, flags }
                    );
                } else {
                    assert!(matches!(result, Err(Error::NonFiniteArgument { .. })));
                }
                assert!(events.is_empty());
            }
        }
    }
}

/// Enforce explicit LEB model-point coverage on the original public epoch.
/// The inclusive fitted interpolated-apsis window precedes stored-channel
/// coverage. Missing model channels are allowed: their analytical models do
/// not require those retired coefficients. The guard runs before configured
/// TT conversion, model evaluation and heliocentric/barycentric zero results.
/// Auto and other provider modes retain their own acquisition policies.
fn model_point_request_guard(
    stack: &ServingStack<'_, '_>,
    ctx: &EphemerisContext,
    jd: f64,
    body: i32,
) -> Result<(), Error> {
    use crate::positions::StateProvider;
    let source_body = crate::lunar::south_node_north_source(body).unwrap_or(body);
    if ctx.mode() != ProviderMode::SealedLeb
        || !matches!(source_body, MEAN_NODE | MEAN_APOG | INTP_APOG | INTP_PERG)
    {
        return Ok(());
    }
    for (start_jd, end_jd) in [
        crate::lunar::interpolated_apse_model_window(source_body),
        stack.serving.provider().body_coverage(source_body),
    ]
    .into_iter()
    .flatten()
    {
        if jd < start_jd || jd > end_jd {
            return Err(Error::EphemerisRange {
                requested_jd: jd,
                start_jd,
                end_jd,
                body: Some(source_body),
            });
        }
    }
    Ok(())
}

/// Live counterpart; callbacks are borrowed only at consumed sidereal stages.
#[allow(clippy::too_many_arguments)]
fn model_point_entry_with_anchors(
    stack: &ServingStack<'_, '_>,
    ctx: &EphemerisContext,
    jd: f64,
    body: i32,
    flags: i32,
    is_ut: bool,
    live: Option<&dyn crate::sidereal::DynamicAnchorProvider>,
) -> Result<([f64; 6], i32), Error> {
    model_point_entry_with_anchors_prepared(stack, ctx, jd, body, flags, is_ut, live, None)
}

#[allow(clippy::too_many_arguments)]
fn model_point_entry_with_anchors_prepared(
    stack: &ServingStack<'_, '_>,
    ctx: &EphemerisContext,
    jd: f64,
    body: i32,
    flags: i32,
    is_ut: bool,
    live: Option<&dyn crate::sidereal::DynamicAnchorProvider>,
    prepared: Option<positions::PreparedCalcFlags>,
) -> Result<([f64; 6], i32), Error> {
    if prepared.is_none() {
        positions::validate_barycentric_moseph(flags, if is_ut { "calc_ut" } else { "calc" })?;
    }
    model_point_request_guard(stack, ctx, jd, body)?;
    match lunar_point_route_with_anchors_prepared(
        stack, ctx, jd, body, flags, is_ut, true, live, prepared,
    ) {
        Err(
            Error::UnknownBody { .. }
            | Error::LebBodyNotFound { .. }
            | Error::LebOutOfRange { .. }
            | Error::EphemerisRange { .. }
            | Error::NonFiniteArgument { .. }
            | Error::SealedChannelNotServed { .. },
        ) => lunar_point_route_with_anchors_prepared(
            stack, ctx, jd, body, flags, is_ut, false, live, prepared,
        ),
        result => result,
    }
}

/// Live counterpart; callbacks are borrowed only at consumed sidereal stages.
#[allow(clippy::too_many_arguments)]
fn true_point_entry_with_anchors(
    stack: &ServingStack<'_, '_>,
    ctx: &EphemerisContext,
    jd: f64,
    body: i32,
    flags: i32,
    is_ut: bool,
    live: Option<&dyn crate::sidereal::DynamicAnchorProvider>,
) -> Result<([f64; 6], i32), Error> {
    true_point_entry_with_anchors_prepared(stack, ctx, jd, body, flags, is_ut, live, None)
}

#[allow(clippy::too_many_arguments)]
fn true_point_entry_with_anchors_prepared(
    stack: &ServingStack<'_, '_>,
    ctx: &EphemerisContext,
    jd: f64,
    body: i32,
    flags: i32,
    is_ut: bool,
    live: Option<&dyn crate::sidereal::DynamicAnchorProvider>,
    prepared: Option<positions::PreparedCalcFlags>,
) -> Result<([f64; 6], i32), Error> {
    use crate::positions::StateProvider;
    let sampled_serving = stack.serving.with_sampled_derivatives();
    let sampled_stack = stack.with_derivative_serving(&sampled_serving);
    let stack = if opaque_live_body_chain(ctx, flags, live) {
        &sampled_stack
    } else {
        stack
    };
    let provider = stack.serving.provider();
    let source_body = crate::lunar::south_node_north_source(body).unwrap_or(body);
    let stored = provider.has_body(source_body);
    let attempted = lunar_point_route_with_anchors_prepared(
        stack, ctx, jd, body, flags, is_ut, stored, live, prepared,
    );
    match attempted {
        Err(
            Error::UnknownBody { .. }
            | Error::LebBodyNotFound { .. }
            | Error::LebOutOfRange { .. }
            | Error::EphemerisRange { .. }
            | Error::NonFiniteArgument { .. }
            | Error::SealedChannelNotServed { .. },
        ) if stored => {
            // Coverage belongs to the original public argument, independently
            // of the configured fast TT that failed. An admitted UT1 request
            // may therefore fall through when its fast TT lies outside K.
            if ctx.mode() == ProviderMode::SealedLeb {
                if let Some((start, end)) = provider.body_coverage(source_body) {
                    if jd < start || jd > end {
                        return Err(Error::EphemerisRange {
                            requested_jd: jd,
                            start_jd: start,
                            end_jd: end,
                            body: Some(source_body),
                        });
                    }
                }
            }
            lunar_point_route_with_anchors_prepared(
                stack, ctx, jd, body, flags, is_ut, false, live, prepared,
            )
        }
        result => result,
    }
}

/// Assemble one fast or native lunar route, retaining that route's TT and
/// sidereal clocks through flag normalization and final output shaping. Native
/// base ayanamsha observations use the independent global UT1→TT conversion;
/// the configured Delta-T provider remains available to defining epochs.
/// Live counterpart; callbacks are borrowed only at consumed sidereal stages.
#[allow(clippy::too_many_arguments)]
fn lunar_point_route_with_anchors_prepared(
    stack: &ServingStack<'_, '_>,
    ctx: &EphemerisContext,
    jd: f64,
    body: i32,
    flags: i32,
    is_ut: bool,
    fast: bool,
    live: Option<&dyn crate::sidereal::DynamicAnchorProvider>,
    prepared: Option<positions::PreparedCalcFlags>,
) -> Result<([f64; 6], i32), Error> {
    let provider = stack.serving.provider();
    let source_body = crate::lunar::south_node_north_source(body).unwrap_or(body);
    let is_model = matches!(source_body, MEAN_NODE | MEAN_APOG | INTP_APOG | INTP_PERG);
    let has_nutation = provider.has_nutation();
    // An absent series selects independent frame data at the same fast TT;
    // an existing series whose sample misses must fail the whole fast attempt.
    let model_nutation: &dyn Fn(f64) -> Result<(f64, f64), Error> = if fast && has_nutation {
        stack.nutation_at
    } else {
        &crate::erfa::nut06a
    };
    let to_tt = |ut| {
        if fast {
            Ok(ut + crate::time::deltat_with_config(ut, stack.delta_cfg)?)
        } else {
            crate::time::default_timescale_tt_from_ut1(ut)
        }
    };
    let raw = |tt, id, calc_flags| {
        if crate::lunar::wants_lunar_helio_bary_zero(id, calc_flags) {
            return Ok(([0.0; 6], calc_flags));
        }
        let values = if matches!(id, MEAN_NODE | MEAN_APOG | INTP_APOG | INTP_PERG) {
            model_point_state(tt, id, calc_flags, fast, model_nutation)?
        } else if fast {
            stored_true_point(stack, ctx, tt, id, calc_flags)?
        } else {
            fallback_true_point(stack, tt, id, calc_flags)?
        };
        Ok((values, calc_flags))
    };
    let run = |compute_tt: &ModelCompute<'_>| {
        let engine = positions::CalcEngine {
            ut1_to_tt_jd: &to_tt,
            compute_tt,
            nutation_obliquity_ut: &positions::ecl_nut_result_ut,
            nutation_obliquity_tt: &|tt, nut_flags| {
                positions::ecl_nut_result(tt, nut_flags, &crate::erfa::nut06a)
            },
            jpl_fallback_active: false,
        };
        let source = if fast {
            "LEB"
        } else if is_model {
            "Analytical"
        } else {
            "Skyfield"
        };
        positions::with_serving_trace_source(body, source, fast, || {
            prepared_engine_entry(&engine, jd, body, flags, is_ut, prepared)
        })
    };
    if flags & FLG_SIDEREAL == 0 {
        return run(&raw);
    }
    let (selection, _) = sidereal_selection_for_context(ctx);
    if !fast
        && resolve_sidereal_route(flags, selection.mode(), selection.bits())
            == Some(SiderealRoute::BaseCorrection)
    {
        // A native fallback samples its base offset at the Time object's UT1.
        // Observation TT is then reconstructed by the independent global
        // timescale, not by the caller's fast-path Delta-T override. Keep the
        // configured provider only for the ayanamsha's defining-epoch inputs.
        let delta = |ut| crate::time::deltat_with_config(ut, stack.delta_cfg);
        let anchors = crate::stars::SealedAnchorProvider::new(provider);
        let anchors = if let Some(live) = live {
            anchors.with_dynamic_override(live)
        } else {
            anchors
        };
        let aya = SealedAyanamshaServing::new(selection, &delta, Some(&anchors));
        let request_tt = if is_ut { to_tt(jd)? } else { jd };
        let sidereal = |tt, id, calc_flags| {
            let (mut value, ret) = raw(tt, id, calc_flags)?;
            if calc_flags & FLG_EQUATORIAL == 0 && value[2] != 0.0 {
                // A Time created from UT1 retains its original UT1 part.
                // Inverting the collapsed scalar TT would quantize that
                // independent clock before the one-second ayanamsha stencil.
                // TT-created side samples instead use their own inverse.
                let ut = if is_ut && tt == request_tt {
                    jd
                } else {
                    crate::time::default_timescale_ut1_from_tt(tt)?
                };
                let offset_at = |ut| {
                    let global_tt = crate::time::default_timescale_tt_from_ut1(ut)?;
                    let (noaberr, nogdefl, use_true) =
                        crate::sidereal::resolve_ayanamsha_flags(calc_flags);
                    if use_true {
                        aya.true_tt_with(global_tt, noaberr, nogdefl)
                    } else {
                        aya.mean_tt_with(global_tt, noaberr, nogdefl)
                    }
                };
                let (offset, _) = offset_at(ut)?;
                let speed = calc_flags & FLG_SPEED != 0;
                let (before, after) = if speed {
                    let half = SIDEREAL_SPEED_HALF_STEP_DAYS;
                    (offset_at(ut - half)?.0, offset_at(ut + half)?.0)
                } else {
                    (offset, offset)
                };
                (value[0], value[3]) =
                    apply_sidereal_correction(value[0], value[3], offset, before, after, speed)?;
            }
            Ok((value, ret))
        };
        return run(&sidereal);
    }
    let delta = |ut| crate::time::deltat_with_config(ut, stack.delta_cfg);
    let delta_value_rate = |tt| {
        ctx.delta_t_value_and_rate_days(tt, stack.delta_cfg.model, stack.delta_cfg.iers_enabled)
    };
    let fast_owner = live.map(|live| FastPositionAnchors {
        live,
        delta: &delta,
        delta_value_rate: &delta_value_rate,
        want_rate: flags & (FLG_SPEED | crate::constants::FLG_SPEED3) != 0,
        clock: std::cell::Cell::new(None),
    });
    let frame_live = if fast {
        fast_owner
            .as_ref()
            .map(|p| p as &dyn crate::sidereal::DynamicAnchorProvider)
    } else {
        live
    };
    // Preserve stored value/derivative ownership unless this route actually
    // replaces the frame feed with the explicit model callback.
    let frame_stack = is_model.then(|| stack.with_nutation(model_nutation));
    serving::with_fixed_epoch_engine_with_anchor_clocks(
        frame_stack.as_ref().unwrap_or(stack),
        ctx,
        &raw,
        frame_live,
        live,
        |engine| run(engine.compute_tt),
    )
}

/// Stored ecliptic-direct state at scalar TT, in degrees/AU and daily rates.
/// Nutation-free side positions use the same channel and stored frame feed.
/// TRUE_NODE's radius follows the osculating conic; missing/degenerate Moon
/// support retains the channel's radius proxy without changing its angles.
fn stored_true_point(
    stack: &ServingStack<'_, '_>,
    ctx: &EphemerisContext,
    tt: f64,
    body: i32,
    flags: i32,
) -> Result<[f64; 6], Error> {
    use crate::positions::StateProvider;
    let provider = stack.serving.provider();
    if provider.serving_coord_type(body) != Some(1) {
        // Custom stored ICRS point channels follow their actual channel type;
        // their component zero must never be interpreted as a longitude.
        return stored_state(stack, ctx, tt, body, flags, false).map(|(value, _)| value);
    }
    let sample = |jd| {
        let (p, _) = provider.eval_raw(body, jd).map_err(|e| {
            positions::escalate_sealed_provider_miss(ctx.mode(), provider, body, jd, e)
        })?;
        if p.len() != 3 {
            return Err(Error::SealedChannelNotServed {
                body_id: body,
                coord_type: provider.serving_coord_type(body).unwrap_or(u32::MAX),
            });
        }
        let mut position = [p[0], p[1], p[2]];
        if body == TRUE_NODE {
            match stored_node_radius(stack, jd) {
                Ok(Some(radius)) => position[2] = radius,
                Ok(None) => {}
                Err(
                    Error::UnknownBody { .. }
                    | Error::LebCorrupt { .. }
                    | Error::LebOutOfRange { .. }
                    | Error::EphemerisRange { .. }
                    | Error::NonFiniteArgument { .. }
                    | Error::SealedChannelNotServed { .. },
                ) => {}
                Err(error) => return Err(error),
            }
        }
        if crate::lunar::no_nutation_frame(flags) {
            let (psi, _) = (stack.nutation_at)(jd)?;
            position[0] = (position[0] - psi.to_degrees()).rem_euclid(360.0);
        }
        Ok(position)
    };
    let central = sample(tt)?;
    let mut value = [central[0], central[1], central[2], 0.0, 0.0, 0.0];
    if flags & FLG_SPEED != 0 {
        let half = if flags & (FLG_SPEED3 | FLG_TOPOCTR) != 0 {
            0.0001
        } else if body == OSCU_APOG {
            0.002
        } else {
            0.05
        };
        match (sample(tt - half), sample(tt + half)) {
            (Ok(before), Ok(after)) => {
                let mut dlon = after[0] - before[0];
                if dlon > 180.0 {
                    dlon -= 360.0;
                } else if dlon < -180.0 {
                    dlon += 360.0;
                }
                value[3] = dlon / (2.0 * half);
                value[4] = (after[1] - before[1]) / (2.0 * half);
                value[5] = (after[2] - before[2]) / (2.0 * half);
            }
            (Err(error), _) | (_, Err(error)) => match error {
                Error::UnknownBody { .. }
                | Error::LebCorrupt { .. }
                | Error::LebOutOfRange { .. }
                | Error::EphemerisRange { .. }
                | Error::NonFiniteArgument { .. }
                | Error::SealedChannelNotServed { .. } => {}
                error => return Err(error),
            },
        }
    }
    lunar_point_frame(value, tt, flags, true, stack.nutation_at)
}

/// Radius in AU of the Moon's instantaneous conic at its ascending node.
/// The scalar stored Moon/Earth states rotate through the fast frame before
/// computing p=h²/GM and e=(v²/GM−1/r)r−(r·v/GM)v. Returns None for a
/// degenerate state, so the direct channel retains its stored radius proxy.
fn stored_node_radius(stack: &ServingStack<'_, '_>, tt: f64) -> Result<Option<f64>, Error> {
    let ((x, y, z), (vx, vy, vz)) =
        crate::lunar::sealed_true_state_at(stack.serving.provider(), stack.nutation_at, tt)?;
    stored_node_radius_from_state(([x, y, z], [vx, vy, vz]))
}

/// Same direct-channel conic reduction after acquired Moon/Earth/frame states.
fn stored_node_radius_from_state(
    ([x, y, z], [vx, vy, vz]): ([f64; 3], [f64; 3]),
) -> Result<Option<f64>, Error> {
    // The direct channel's conic uses the Gaussian/Sun-to-Earth-Moon
    // parameter, distinct from the native lunar geometry's DE430 GM.
    let gm = (0.01720209895_f64 * 0.01720209895_f64) / 328900.5596;
    let hx = y * vz - z * vy;
    let hy = z * vx - x * vz;
    let hz = x * vy - y * vx;
    let p = (hx * hx + hy * hy + hz * hz) / gm;
    let r = (x * x + y * y + z * z).sqrt();
    if r == 0.0 {
        return Ok(None);
    }
    let c1 = (vx * vx + vy * vy + vz * vz) / gm - 1.0 / r;
    let c2 = (x * vx + y * vy + z * vz) / gm;
    let ex = c1 * x - c2 * vx;
    let ey = c1 * y - c2 * vy;
    let ez = c1 * z - c2 * vz;
    let e = (ex * ex + ey * ey + ez * ez).sqrt();
    let n = (-hy).hypot(hx);
    if n <= 0.0 || e <= 0.0 {
        return Ok(None);
    }
    let cosine = ((-hy * ex + hx * ey) / (n * e)).clamp(-1.0, 1.0);
    let denominator = 1.0 + e * cosine;
    if denominator == 0.0 {
        return Ok(None);
    }
    Ok(Some(p / denominator))
}

/// Missing-channel lunar fallback over scalar LEB Moon/Earth states. The
/// native Time frame uses P03 precession at TDB and IAU 2000A nutation at TT;
/// neither caller Delta-T nor stored fast-frame nutation replaces those inputs.
fn fallback_true_point(
    stack: &ServingStack<'_, '_>,
    tt: f64,
    body: i32,
    flags: i32,
) -> Result<[f64; 6], Error> {
    use crate::positions::{StateProvider, reference_frames as sf};
    let provider = stack.serving.provider();
    // Resolve both endpoint aliases before evaluating either state. The
    // Earth alias owns a missing-endpoint refusal even when Moon is missing
    // or its requested sample lies outside coverage; vector reads stay Moon
    // first after successful alias resolution.
    for endpoint in [EARTH, MOON] {
        if !provider.has_body(endpoint) {
            return Err(Error::UnknownBody { body_id: endpoint });
        }
    }
    let state_at = |jd: f64| {
        let relative = provider
            .relative_geometric_state(MOON, EARTH, jd)
            .map_err(|error| {
                let failed_body = match error.primary_error() {
                    Error::LebBodyNotFound { body_id } | Error::LebOutOfRange { body_id, .. } => {
                        *body_id
                    }
                    _ => MOON,
                };
                positions::escalate_sealed_provider_miss(
                    ProviderMode::SealedLeb,
                    provider,
                    failed_body,
                    jd,
                    error,
                )
            })?;
        let whole = jd.floor();
        let fraction = jd - whole;
        let tdb = whole + (fraction + sf::tdb_minus_tt_seconds(whole, fraction)? / 86400.0);
        let eps = sf::mean_obliquity_arcsec(tdb)? * sf::FRAME_ASEC2RAD;
        let (psi, deps) = sf::iau2000a_tenth_microarcsec(jd)?;
        let psi = psi * sf::FRAME_TENTH_USEC_2_RAD;
        let deps = deps * sf::FRAME_TENTH_USEC_2_RAD;
        let nut = sf::nutation_matrix(eps, eps + deps, psi)?;
        let equinox = sf::equinox_matrix(nut, sf::precession_matrix(tdb)?)?;
        let frame = sf::ecliptic_date_matrix(equinox, eps, deps)?;
        let r = sf::matvec3(frame, relative.pos())?;
        let v = sf::matvec3(frame, relative.vel())?;
        Ok(((r[0], r[1], r[2]), (v[0], v[1], v[2])))
    };
    let nutation = crate::erfa::nut06a;
    let dpsi = |jd| nutation(jd).map(|(psi, _)| psi);
    let rate =
        |jd: f64, half: f64| Ok((dpsi(jd + half)? - dpsi(jd - half)?).to_degrees() / (2.0 * half));
    let value =
        crate::lunar::true_lunar_state_with_providers(body, tt, flags, &state_at, &dpsi, &rate)?;
    lunar_point_frame(value, tt, flags, false, &nutation)
}

/// Map ecliptic degrees/AU and daily rates into the requested epoch/equator.
/// Fast direct channels and models use the 0.001-day centered mapping; native
/// fallback J2000 precession uses 0.5 day and equatorial rates use cotrans_sp.
/// Fast requests without SPEED discard mapped rates before output shaping.
/// Radius and radial rate are invariant. Frame acquisition failures propagate.
fn lunar_point_frame(
    mut value: [f64; 6],
    tt: f64,
    flags: i32,
    fast: bool,
    nutation_at: &dyn Fn(f64) -> Result<(f64, f64), Error>,
) -> Result<[f64; 6], Error> {
    if flags & (FLG_J2000 | FLG_EQUATORIAL) == 0 {
        return Ok(value);
    }
    let eps = if flags & FLG_J2000 != 0 {
        84381.406 / 3600.0
    } else {
        positions::vondrak_mean_obliquity_deg(tt)?
            + if flags & (FLG_NONUT | FLG_SIDEREAL) == 0 {
                nutation_at(tt)?.1.to_degrees()
            } else {
                0.0
            }
    };
    let mapped = |lon, lat, epoch| {
        let (lon, lat) = if flags & FLG_J2000 != 0 {
            crate::hypothetical::precess_ecliptic(lon, lat, epoch, 2451545.0)?
        } else {
            (lon, lat)
        };
        if flags & FLG_EQUATORIAL != 0 && fast {
            positions::cotrans(lon, lat, -eps)
        } else {
            Ok((lon, lat))
        }
    };
    let (lon, lat) = mapped(value[0], value[1], tt)?;
    if fast || (flags & FLG_J2000 != 0 && flags & FLG_SPEED != 0) {
        let half = if fast { 0.001 } else { 0.5 };
        let (lm, bm) = mapped(
            value[0] - value[3] * half,
            value[1] - value[4] * half,
            tt - half,
        )?;
        let (lp, bp) = mapped(
            value[0] + value[3] * half,
            value[1] + value[4] * half,
            tt + half,
        )?;
        let mut dl = lp - lm;
        if dl > 180.0 {
            dl -= 360.0;
        } else if dl < -180.0 {
            dl += 360.0;
        }
        value[3] = dl / (2.0 * half);
        value[4] = (bp - bm) / (2.0 * half);
    }
    value[0] = lon;
    value[1] = lat;
    if !fast && flags & FLG_EQUATORIAL != 0 {
        let (p, v) = crate::utils::cotrans_sp(
            value[0], value[1], value[2], value[3], value[4], value[5], -eps,
        )?;
        value = [p[0], p[1], p[2], v[0], v[1], v[2]];
    }
    if fast && flags & FLG_SPEED == 0 {
        value[3..].fill(0.0);
    }
    Ok(value)
}

/// Select the prepared or ordinary engine entry without reinterpreting flags.
/// The original epoch scale chooses TT or UT1 in either case; clock conversion,
/// normalization and numerical errors remain owned by that selected entry.
fn prepared_engine_entry(
    engine: &positions::CalcEngine<'_>,
    jd: f64,
    body: i32,
    flags: i32,
    is_ut: bool,
    prepared: Option<positions::PreparedCalcFlags>,
) -> Result<([f64; 6], i32), Error> {
    if let Some(prepared) = prepared {
        if is_ut {
            positions::calc_ut_prepared(engine, jd, body, prepared)
        } else {
            positions::calc_prepared(engine, jd, body, prepared)
        }
    } else if is_ut {
        positions::calc_ut(engine, jd, body, flags)
    } else {
        positions::calc(engine, jd, body, flags)
    }
}

#[cfg(test)]
fn curated_entry(
    stack: &ServingStack<'_, '_>,
    ctx: &EphemerisContext,
    jd: f64,
    body: i32,
    flags: i32,
    is_ut: bool,
) -> Result<([f64; 6], i32), Error> {
    curated_entry_with_anchors(stack, ctx, jd, body, flags, is_ut, None)
}

/// Curated-minor serving over one borrowed sealed support stack.
/// The complete stored-state attempt uses configured fast TT, stored frames
/// and direct Cartesian output. A recoverable absent/range/frame-channel miss
/// retries the complete local model at independent global TT. The companion
/// body's stored range never bounds that fallback; Earth/Sun support still
/// does. [`Error::EphemerisRange`] is terminal and retains the requested body;
/// raw channel misses can still retry when the support covers the epoch.
/// Corrupt channels and missing-observer preconditions remain fatal.
/// Explicit LEB neither reads JPL nor applies the non-LEB strict-SPK policy.
/// The UT stored branch retains the caller's central UT1 and freezes its
/// TT−UT1 difference for side samples; the TT branch obtains UT1 with the
/// fast path's separate one-step fixed-point inverse. Live callbacks are
/// borrowed only at consumed sidereal stages.
#[allow(clippy::too_many_arguments)]
fn curated_entry_with_anchors(
    stack: &ServingStack<'_, '_>,
    ctx: &EphemerisContext,
    jd: f64,
    body: i32,
    flags: i32,
    is_ut: bool,
    live: Option<&dyn crate::sidereal::DynamicAnchorProvider>,
) -> Result<([f64; 6], i32), Error> {
    curated_entry_with_anchors_prepared(stack, ctx, jd, body, flags, is_ut, live, None)
}

#[allow(clippy::too_many_arguments)]
fn curated_entry_with_anchors_prepared(
    stack: &ServingStack<'_, '_>,
    ctx: &EphemerisContext,
    jd: f64,
    body: i32,
    flags: i32,
    is_ut: bool,
    live: Option<&dyn crate::sidereal::DynamicAnchorProvider>,
    prepared: Option<positions::PreparedCalcFlags>,
) -> Result<([f64; 6], i32), Error> {
    curated_entry_core_prepared(stack, ctx, jd, body, flags, is_ut, live, true, prepared)
}

/// Shared curated attempt with either the original typed or actual-error catch.
#[allow(clippy::too_many_arguments)]
fn curated_entry_core_prepared(
    stack: &ServingStack<'_, '_>,
    ctx: &EphemerisContext,
    jd: f64,
    body: i32,
    flags: i32,
    is_ut: bool,
    live: Option<&dyn crate::sidereal::DynamicAnchorProvider>,
    retry_inside: bool,
    prepared: Option<positions::PreparedCalcFlags>,
) -> Result<([f64; 6], i32), Error> {
    let sampled_serving = stack.serving.with_sampled_derivatives();
    let sampled_stack = stack.with_derivative_serving(&sampled_serving);
    let stack = if opaque_live_body_chain(ctx, flags, live) {
        &sampled_stack
    } else {
        stack
    };
    // Raw flag validation precedes the fast observer precondition; observer
    // configuration in turn precedes configured time conversion and data.
    if prepared.is_none() {
        positions::validate_barycentric_moseph(flags, if is_ut { "calc_ut" } else { "calc" })?;
    }
    let site = ctx
        .observer()
        .map(|o| (o.lon_deg(), o.lat_deg(), o.height_m()));
    positions::resolve_topo_request(
        prepared.map_or_else(
            || positions::normalize_calc_flags(flags, false),
            |p| p.normalized,
        ),
        positions::remap_ast_offset(body),
        site,
        None,
    )?;
    let delta = |ut| crate::time::deltat_with_config(ut, stack.delta_cfg);
    let site_ut1 = |tt, calc_flags| {
        if calc_flags & FLG_TOPOCTR != 0 {
            fast_ut1_partner(tt, jd, is_ut, &delta).map(Some)
        } else {
            Ok(None)
        }
    };
    let raw = |tt, id, calc_flags| {
        stored_state_with_ut1(
            stack,
            ctx,
            tt,
            id,
            calc_flags,
            false,
            site_ut1(tt, calc_flags)?,
        )
    };
    let to_tt = |ut| Ok(ut + crate::time::deltat_with_config(ut, stack.delta_cfg)?);
    let delta_value_rate = |tt| {
        ctx.delta_t_value_and_rate_days(tt, stack.delta_cfg.model, stack.delta_cfg.iers_enabled)
    };
    let fast_owner = live.map(|live| FastPositionAnchors {
        live,
        delta: &delta,
        delta_value_rate: &delta_value_rate,
        want_rate: flags & (FLG_SPEED | crate::constants::FLG_SPEED3) != 0,
        clock: std::cell::Cell::new(None),
    });
    let fast_live = fast_owner
        .as_ref()
        .map(|p| p as &dyn crate::sidereal::DynamicAnchorProvider);
    let attempted = positions::with_serving_trace_source(body, "LEB", true, || {
        if flags & FLG_SIDEREAL != 0 {
            serving::with_fixed_epoch_engine_with_anchor_clocks(
                stack,
                ctx,
                &raw,
                fast_live,
                live,
                |engine| prepared_engine_entry(engine, jd, body, flags, is_ut, prepared),
            )
        } else {
            let engine = positions::CalcEngine {
                ut1_to_tt_jd: &to_tt,
                compute_tt: &raw,
                nutation_obliquity_ut: &positions::ecl_nut_result_ut,
                nutation_obliquity_tt: &|tt, nut_flags| {
                    positions::ecl_nut_result(tt, nut_flags, &crate::erfa::nut06a)
                },
                jpl_fallback_active: false,
            };
            let cartesian = |tt, id, calc_flags| {
                // The observer precondition must precede the absent-body signal
                // even when XYZ selects the direct Cartesian fast arm.
                match stored_state_with_ut1(
                    stack,
                    ctx,
                    tt,
                    id,
                    calc_flags,
                    true,
                    site_ut1(tt, calc_flags)?,
                ) {
                    Err(Error::SealedChannelNotServed { .. }) => Ok(None),
                    result => result.map(Some),
                }
            };
            if let Some(prepared) = prepared {
                positions::calc_with_cartesian_prepared(
                    &engine, &cartesian, jd, body, prepared, is_ut,
                )
            } else if is_ut {
                positions::calc_ut_with_cartesian(&engine, &cartesian, jd, body, flags)
            } else {
                positions::calc_with_cartesian(&engine, &cartesian, jd, body, flags)
            }
        }
    });
    if !retry_inside {
        return attempted;
    }
    match attempted {
        Err(
            Error::UnknownBody { .. }
            | Error::LebBodyNotFound { .. }
            | Error::LebOutOfRange { .. }
            | Error::NonFiniteArgument { .. }
            | Error::SealedChannelNotServed { .. },
        ) => model_entry_with_models_and_anchors(None, stack, ctx, jd, body, flags, is_ut, live),
        result => result,
    }
}

/// UT1 partner of a stored-state TT sample, in Julian days.
/// An original UT request retains that exact value at its configured central
/// TT; additional callback samples subtract the frozen central TT−UT1 pair.
/// A TT request uses the canonical fast inverse with one fixed-point step,
/// independently of the global Time object used by the Keplerian retry.
/// The configured Delta-T callback consumes UT1 and returns days; its errors
/// propagate. Provider-side SPEED stencils freeze this returned central pair.
fn fast_ut1_partner(
    tt: f64,
    request_jd: f64,
    is_ut: bool,
    delta_t: &dyn Fn(f64) -> Result<f64, Error>,
) -> Result<f64, Error> {
    if is_ut {
        let central_tt = request_jd + delta_t(request_jd)?;
        if tt == central_tt {
            Ok(request_jd)
        } else {
            Ok(tt - (central_tt - request_jd))
        }
    } else {
        let initial_delta = delta_t(tt)?;
        Ok(tt - delta_t(tt - initial_delta)?)
    }
}

/// Tropical dispatch at scalar TT over stored data. Ordinary bodies use the
/// existing topocentric/geocentric operations. The lazy clock is consumed
/// only by an actual site correction, after observer preconditions. Model points consume stored
/// nutation and the fast ecliptic frame mapping; true points have a separate
/// direct/fallback acquisition dispatcher.
fn tropical(
    stack: &ServingStack<'_, '_>,
    ctx: &EphemerisContext,
    tt: f64,
    body: i32,
    flags: i32,
    site_clock: &SiteClock<'_>,
    public_epoch: f64,
) -> Result<([f64; 6], i32), Error> {
    if !crate::lunar::is_lunar_point(body) {
        return stored_state_with_clock_at(
            stack,
            ctx,
            tt,
            body,
            flags,
            false,
            site_clock,
            Some(public_epoch),
        );
    }
    if crate::lunar::wants_lunar_helio_bary_zero(body, flags) {
        return Ok(([0.0; 6], flags));
    }
    model_point_state(tt, body, flags, true, stack.nutation_at).map(|values| (values, flags))
}

/// Analytical mean/interpolated state at TT in degrees/AU and daily rates.
/// The computation flags retain SIDEREAL/EQUATORIAL for the mean-ecliptic
/// nutation rule; epoch/equator mapping follows state acquisition, and the
/// public CalcEngine owns representation shaping and south-node recursion.
/// Fast frames use the 0.001-day angular mapping with stored nutation; native
/// frames use the independent nutation feed and native precession/rate tail.
/// Fast clean models retain their nutation-rate acquisition without SPEED
/// before discarding rates: mean true-ecliptic and interpolated mean-ecliptic
/// requests still require both 0.001-day side frame samples. Native requests
/// acquire this rate only under SPEED. Central frame acquisition comes first.
fn model_point_state(
    tt: f64,
    body: i32,
    flags: i32,
    fast: bool,
    nutation_at: &dyn Fn(f64) -> Result<(f64, f64), Error>,
) -> Result<[f64; 6], Error> {
    let dpsi = |jd| nutation_at(jd).map(|(dpsi, _)| dpsi);
    let rate = |jd: f64, half: f64| {
        if !half.is_finite() || half <= 0.0 {
            return Err(Error::NonFiniteArgument { value: half });
        }
        let before = dpsi(jd - half)?;
        let after = dpsi(jd + half)?;
        Ok((after - before).to_degrees() / (2.0 * half))
    };
    let values = match body {
        MEAN_NODE | MEAN_APOG => {
            let (a, b, c, d, e, f) = crate::lunar::mean_lunar_point_state(
                body,
                tt,
                flags,
                crate::lunar::LunarCoverage::Unenforced,
                &dpsi,
                &rate,
            )?;
            [a, b, c, d, e, f]
        }
        INTP_APOG | INTP_PERG => {
            let (a, b, c, d, e, f) =
                crate::lunar::interp_lunar_point_state(body, tt, flags, &dpsi, &rate)?;
            [a, b, c, d, e, f]
        }
        _ => return Err(Error::UnknownBody { body_id: body }),
    };
    if fast && flags & FLG_SPEED == 0 {
        let no_nutation = crate::lunar::no_nutation_frame(flags);
        if (crate::lunar::is_mean_lunar_point(body) && !no_nutation)
            || (crate::lunar::is_interp_lunar_point(body) && no_nutation)
        {
            // The fast pipeline reads this rate even though its later
            // no-SPEED tail discards it. A side miss must retry the whole
            // route rather than retaining configured fast TT and its frame.
            rate(tt, 0.001)?;
        }
    }
    lunar_point_frame(values, tt, flags, fast, nutation_at)
}

/// Direct Cartesian provider, with the same observer/UT partner as spherical
/// serving. Preserves source misses and all native correction/frame operations.
fn stored_state(
    stack: &ServingStack<'_, '_>,
    ctx: &EphemerisContext,
    tt: f64,
    body: i32,
    flags: i32,
    xyz: bool,
) -> Result<([f64; 6], i32), Error> {
    stored_state_with_ut1(stack, ctx, tt, body, flags, xyz, None)
}

/// Stored-state computation with an optional caller-owned UT1 site partner.
/// `Some` preserves the public entry's UT1/TT pairing and the topocentric
/// core's frozen side-sample policy. `None` retains the existing configured
/// inverse used by other callers. Spherical/Cartesian provider operations,
/// observer and coverage errors retain the [`stored_state`] contract.
fn stored_state_with_ut1(
    stack: &ServingStack<'_, '_>,
    ctx: &EphemerisContext,
    tt: f64,
    body: i32,
    flags: i32,
    xyz: bool,
    site_ut1: Option<f64>,
) -> Result<([f64; 6], i32), Error> {
    let site_clock = |epoch| match site_ut1 {
        Some(ut) => Ok(ut),
        None => Ok(epoch - crate::time::deltat_with_config(epoch, stack.delta_cfg)?),
    };
    stored_state_with_clock(stack, ctx, tt, body, flags, xyz, &site_clock)
}

/// Stored spherical or Cartesian state with a lazy UT1 site clock.
/// Observer resolution precedes clock evaluation. Geocentric and lunar
/// no-parallax centers never call the clock. A topocentric center evaluates
/// the clock once at its central TT, then the provider freezes TT−UT1 for
/// SPEED side samples. Clock errors propagate before provider acquisition;
/// stored coverage, frames, WGS84 geometry and public shaping remain shared.
fn stored_state_with_clock(
    stack: &ServingStack<'_, '_>,
    ctx: &EphemerisContext,
    tt: f64,
    body: i32,
    flags: i32,
    xyz: bool,
    site_clock: &SiteClock<'_>,
) -> Result<([f64; 6], i32), Error> {
    stored_state_with_clock_at(stack, ctx, tt, body, flags, xyz, site_clock, None)
}

/// Same acquisition with the original public TT/UT1 retained for a consumed
/// sealed reader refusal. The provider failure is captured before semantic
/// escalation; metadata comes from that same owned serving view.
#[allow(clippy::too_many_arguments)]
fn stored_state_with_clock_at(
    stack: &ServingStack<'_, '_>,
    ctx: &EphemerisContext,
    tt: f64,
    body: i32,
    flags: i32,
    xyz: bool,
    site_clock: &SiteClock<'_>,
    public_epoch: Option<f64>,
) -> Result<([f64; 6], i32), Error> {
    let site = ctx
        .observer()
        .map(|o| (o.lon_deg(), o.lat_deg(), o.height_m()));
    let result = match positions::resolve_topo_request(flags, body, site, None)? {
        TopoPlan::Geocentric | TopoPlan::TopoNoop => {
            if xyz {
                positions::serving_compute_xyz_tt(stack.serving, tt, body, flags)
            } else {
                positions::serving_compute_tt(stack.serving, tt, body, flags)
            }
        }
        TopoPlan::Topocentric {
            lon_deg,
            lat_deg,
            height_m,
        } => {
            let ut = site_clock(tt)?;
            let rotation = |t, u| crate::erfa::c2t06a(t, u, 0.0, 0.0);
            if xyz {
                positions::serving_compute_xyz_tt_topocentric(
                    stack.serving,
                    tt,
                    ut,
                    body,
                    flags,
                    (lon_deg, lat_deg, height_m),
                    &rotation,
                    // Scalar `fast_calc` behavior: no surface-observer
                    // Earth-mass term (the phenomena vector branch opts
                    // in separately at its own call site).
                    false,
                    // Scalar `fast_calc` site model: WGS84 (the vector
                    // backend's IERS2010 `Topos` model stays on the
                    // phenomena vector arm).
                    positions::SiteEllipsoid::Wgs84,
                )
            } else {
                positions::serving_compute_tt_topocentric(
                    stack.serving,
                    tt,
                    ut,
                    body,
                    flags,
                    (lon_deg, lat_deg, height_m),
                    &rotation,
                    // Scalar `fast_calc` behavior: no surface-observer
                    // Earth-mass term (the phenomena vector branch opts
                    // in separately at its own call site).
                    false,
                    // Scalar `fast_calc` site model: WGS84 (the vector
                    // backend's IERS2010 `Topos` model stays on the
                    // phenomena vector arm).
                    positions::SiteEllipsoid::Wgs84,
                )
            }
        }
    };
    // Both observer routes expose requested-body coverage errors, not raw
    // support-body misses. Auto retains its existing non-sealed taxonomy.
    result.map_err(|e| {
        let provider = stack.serving.provider();
        let public_error = public_epoch.and_then(|epoch| {
            (ctx.mode() == ProviderMode::SealedLeb)
                .then(|| provider.public_calc_range_error(&e, body, epoch))
                .flatten()
        });
        match public_error {
            Some(error) => error,
            None => positions::escalate_sealed_provider_miss(ctx.mode(), provider, body, tt, e),
        }
    })
}

#[cfg(test)]
mod curated_clock_tests {
    use super::*;
    use std::cell::RefCell;

    #[test]
    fn observer_and_raw_center_errors_precede_unavailable_delta_t() {
        use sha2::{Digest, Sha256};

        let path = std::env::var_os("LIBEPHEMERIS_TEST_LEB2")
            .map(std::path::PathBuf::from)
            .unwrap_or_else(crate::test_data::core_leb_path);
        let bytes = std::fs::read(&path).expect("provision the pinned historical fixture");
        let digest: String = Sha256::digest(&bytes)
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect();
        assert!(crate::test_data::admitted_core_hash(&(digest)));
        let mut ctx = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        ctx.attach_leb_file(
            crate::leb::LebFile::from_bytes(path.display().to_string(), bytes).unwrap(),
        );
        serving::with_sealed_serving_stack(
            &ctx,
            crate::time::DeltaTModel::Smh2016,
            false,
            true,
            |stack| {
                let unavailable = crate::time::DeltaTConfig::default();
                let no_frames = |_: f64| -> Result<(f64, f64), Error> {
                    panic!("preflight errors precede stored frame acquisition")
                };
                let mut missing_time = stack.with_nutation(&no_frames);
                missing_time.delta_cfg = &unavailable;
                for (is_ut, entry) in [(true, "calc_ut"), (false, "calc")] {
                    assert!(matches!(
                        curated_entry(&missing_time, &ctx, 2451545.125, CERES, FLG_TOPOCTR, is_ut),
                        Err(Error::ConfigurationMissing {
                            missing_config: "observer_location",
                            ..
                        })
                    ));
                    let flags = FLG_BARYCTR | FLG_MOSEPH;
                    assert_eq!(
                        curated_entry(&missing_time, &ctx, 2451545.125, CERES, flags, is_ut),
                        Err(Error::UnsupportedFlagCombination { entry, flags })
                    );
                    // TOPOCTR is the selected center when it conflicts with
                    // BARYCTR, so its missing observer retains precedence.
                    assert!(matches!(
                        curated_entry(
                            &missing_time,
                            &ctx,
                            2451545.125,
                            CERES,
                            flags | FLG_TOPOCTR,
                            is_ut
                        ),
                        Err(Error::ConfigurationMissing {
                            missing_config: "observer_location",
                            ..
                        })
                    ));
                }
                Ok(())
            },
        )
        .unwrap();
    }

    #[test]
    fn ut_site_partner_preserves_caller_word_and_frozen_side_difference() {
        let ut = 2451545.125;
        let delta = |jd: f64| Ok(0.25 + (jd - ut) / 8.0);
        let tt = ut + 0.25;
        let central = fast_ut1_partner(tt, ut, true, &delta).unwrap();
        assert_eq!(central.to_bits(), ut.to_bits());
        for sample in [tt - 0.001, tt + 0.001] {
            let partner = fast_ut1_partner(sample, ut, true, &delta).unwrap();
            assert_eq!(partner.to_bits(), (sample - 0.25).to_bits());
        }
        // An inverse evaluation at TT would move this deliberately varying
        // clock and would therefore change the central rotation epoch.
        assert_ne!(central, tt - delta(tt).unwrap());
    }

    #[test]
    fn tt_site_partner_uses_one_fixed_point_step_in_configured_delta_t() {
        let tt = 2451545.125;
        let calls = RefCell::new(vec![]);
        let delta = |jd: f64| {
            calls.borrow_mut().push(jd);
            Ok(0.25 + (jd - tt) / 8.0)
        };
        let partner = fast_ut1_partner(tt, tt, false, &delta).unwrap();
        assert_eq!(&*calls.borrow(), &[tt, tt - 0.25]);
        assert_eq!(partner.to_bits(), (tt - 0.21875).to_bits());
    }
}

#[cfg(test)]
mod leb_site_clock_tests;

#[cfg(test)]
mod leb_canonical_tt_fallback_tests;

#[allow(clippy::too_many_arguments)]
fn model_entry_with_models_and_anchors(
    models: Option<&positions::AnalyticalModelCallback<'_>>,
    stack: &ServingStack<'_, '_>,
    ctx: &EphemerisContext,
    jd: f64,
    body: i32,
    flags: i32,
    is_ut: bool,
    live: Option<&dyn crate::sidereal::DynamicAnchorProvider>,
) -> Result<([f64; 6], i32), Error> {
    model_entry_with_models_and_anchors_prepared(
        models, stack, ctx, jd, body, flags, is_ut, live, None,
    )
}

#[allow(clippy::too_many_arguments)]
fn model_entry_with_models_and_anchors_prepared(
    models: Option<&positions::AnalyticalModelCallback<'_>>,
    stack: &ServingStack<'_, '_>,
    ctx: &EphemerisContext,
    jd: f64,
    body: i32,
    flags: i32,
    is_ut: bool,
    live: Option<&dyn crate::sidereal::DynamicAnchorProvider>,
    prepared: Option<positions::PreparedCalcFlags>,
) -> Result<([f64; 6], i32), Error> {
    let (selection, _) = sidereal_selection_for_context(ctx);
    let delta = |ut| crate::time::deltat_with_config(ut, stack.delta_cfg);
    let anchors = crate::stars::SealedAnchorProvider::new(stack.serving.provider());
    let anchors = if let Some(live) = live {
        anchors.with_dynamic_override(live)
    } else {
        anchors
    };
    let aya = SealedAyanamshaServing::new(selection, &delta, Some(&anchors));
    let curated =
        crate::minor_bodies::minor_body_elements(positions::remap_ast_offset(body)).is_some();
    let offset = |ut, calc_flags| {
        if curated || live.is_some() {
            // Local curated output samples its Time object's UT1, then the
            // global timescale supplies observation TT for the ayanamsha.
            // Configured Delta-T remains available to defining epochs.
            let tt = crate::time::default_timescale_tt_from_ut1(ut)?;
            let (noaberr, nogdefl, use_true) = crate::sidereal::resolve_ayanamsha_flags(calc_flags);
            let (value, _) = if use_true {
                aya.true_tt_with(tt, noaberr, nogdefl)?
            } else {
                aya.mean_tt_with(tt, noaberr, nogdefl)?
            };
            Ok(value)
        } else {
            aya.for_flags_ut(ut, calc_flags).map(|(value, _)| value)
        }
    };
    let site = ctx
        .observer()
        .map(|o| (o.lon_deg(), o.lat_deg(), o.height_m()));
    let raw = |global_tt: f64, id, calc_flags| {
        let request_ut1 = if is_ut {
            jd
        } else {
            crate::time::default_timescale_ut1_from_tt(global_tt)?
        };
        positions::context_serving::model_position_with_owned_support(
            models,
            stack.serving.provider(),
            request_ut1,
            global_tt,
            id,
            calc_flags,
            site,
            &delta,
            &offset,
            selection.mode(),
        )
    };
    let to_tt = crate::time::default_timescale_tt_from_ut1;
    let run = |compute_tt: &ModelCompute<'_>| {
        let engine = positions::CalcEngine {
            ut1_to_tt_jd: &to_tt,
            compute_tt,
            nutation_obliquity_ut: &positions::ecl_nut_result_ut,
            nutation_obliquity_tt: &|tt, nut_flags| {
                positions::ecl_nut_result(tt, nut_flags, &crate::erfa::nut06a)
            },
            jpl_fallback_active: false,
        };
        let source = if curated { "Keplerian" } else { "Analytical" };
        positions::with_serving_trace_source(body, source, false, || {
            prepared_engine_entry(&engine, jd, body, flags, is_ut, prepared)
        })
    };
    let route = resolve_sidereal_route(flags, selection.mode(), selection.bits());
    if matches!(route, None | Some(SiderealRoute::BaseCorrection)) {
        run(&raw)
    } else {
        serving::with_fixed_epoch_engine_with_anchors(stack, ctx, &raw, live, |engine| {
            run(engine.compute_tt)
        })
    }
}

/// Whole-source eligibility for live base-correction derivatives. Fixed-frame
/// projections and formula modes do not consume a live anchor and stay separate.
fn opaque_live_body_chain(
    ctx: &EphemerisContext,
    flags: i32,
    live: Option<&dyn crate::sidereal::DynamicAnchorProvider>,
) -> bool {
    let Some(live) = live else { return false };
    let (selection, _) = sidereal_selection_for_context(ctx);
    let mode = strip_to_base_mode(selection.mode());
    crate::sidereal::DYNAMIC_MODE_ANCHORS
        .iter()
        .any(|(id, key)| *id == mode && crate::sidereal::anchor_star(key).is_some())
        && resolve_sidereal_route(flags, mode, selection.bits())
            == Some(crate::sidereal::SiderealRoute::BaseCorrection)
        && !live.has_analytic_anchor_rates()
}

#[cfg(test)]
mod opaque_anchor_chain_tests {
    use super::*;
    use crate::stars::{
        AcquiredAnchorRow, AnchorRowAcquisition, AnchorSample, AnchorSampleOptions,
        NativeStarInputs,
    };
    use std::sync::{Arc, Mutex};
    struct Row {
        times: Arc<Mutex<Vec<f64>>>,
        refuse: bool,
    }
    fn marker() -> Error {
        Error::MathDomain { value: 17.0 }
    }
    fn longitude(tt: f64) -> f64 {
        200.0 + (tt - 2451545.0) * 0.01 + (tt - 2451545.0).powi(2) * 1e-5
    }
    impl AcquiredAnchorRow for Row {
        fn sample(&self, tt: f64, _: AnchorSampleOptions) -> Result<AnchorSample, Error> {
            self.times.lock().unwrap().push(tt);
            if self.refuse {
                Err(marker())
            } else {
                Ok(AnchorSample::Cached(longitude(tt)))
            }
        }
        fn galactic_inputs(&self) -> Result<NativeStarInputs, Error> {
            panic!("Spica must not acquire galactic constructor inputs")
        }
    }
    struct Rows {
        times: Arc<Mutex<Vec<f64>>>,
        refuse: bool,
    }
    impl AnchorRowAcquisition for Rows {
        fn acquire_row(&self, key: &str) -> Result<Arc<dyn AcquiredAnchorRow>, Error> {
            assert_eq!(key, "SPICA");
            Ok(Arc::new(Row {
                times: self.times.clone(),
                refuse: self.refuse,
            }))
        }
    }
    #[test]
    fn opaque_anchor_raw_leb_keeps_complete_sampled_body_and_callback_order() {
        use sha2::{Digest, Sha256};
        let path = std::env::var_os("LIBEPHEMERIS_TEST_LEB2")
            .map(std::path::PathBuf::from)
            .unwrap_or_else(crate::test_data::core_leb_path);
        let bytes = std::fs::read(&path).unwrap();
        assert!(crate::test_data::admitted_core_hash(
            &(Sha256::digest(&bytes)
                .iter()
                .map(|byte| format!("{byte:02x}"))
                .collect::<String>())
        ));
        let mut ctx = EphemerisContext::new()
            .with_mode(ProviderMode::SealedLeb)
            .with_user_delta_t(Some(123.0 / 86400.0));
        ctx.set_sidereal_mode(crate::constants::SIDM_TRUE_CITRA, 0.0, 0.0);
        ctx.attach_leb_file(
            crate::leb::LebFile::from_bytes(path.display().to_string(), bytes).unwrap(),
        );
        let tt = 2451545.125;
        let flags = FLG_SPEED | FLG_SIDEREAL | crate::constants::FLG_SWIEPH;
        let times = Arc::new(Mutex::new(Vec::new()));
        let rows = Rows {
            times: times.clone(),
            refuse: false,
        };
        let session = Session::from_context(ctx.clone());
        let actual = session
            .calc_with_acquisitions(
                tt,
                MARS,
                flags,
                &PositionAcquisitions {
                    anchor_rows: Some(&rows),
                    ..Default::default()
                },
            )
            .unwrap();
        let expected = serving::with_sealed_serving_stack(
            &ctx,
            crate::time::DeltaTModel::Smh2016,
            false,
            true,
            |stack| {
                let sampled = stack.serving.with_sampled_derivatives();
                let (mut xx, _) =
                    positions::serving_compute_tt(&sampled, tt, MARS, flags & !FLG_SIDEREAL)?;
                let at = |at: f64| -> Result<f64, Error> {
                    let g = crate::time::default_timescale_tt_from_ut1(at - 123.0 / 86400.0)?;
                    Ok(longitude(g) - 180.0 + (stack.nutation_at)(at)?.0.to_degrees())
                };
                xx[0] = crate::sidereal::apply_sidereal_offset_deg(xx[0], at(tt)?)?;
                let h = crate::sidereal::SIDEREAL_SPEED_HALF_STEP_DAYS;
                xx[3] -=
                    crate::sidereal::sidereal_rate_correction_deg_day(at(tt - h)?, at(tt + h)?)?;
                Ok(xx)
            },
        )
        .unwrap();
        assert_eq!(actual.values.map(f64::to_bits), expected.map(f64::to_bits));
        let h = crate::sidereal::SIDEREAL_SPEED_HALF_STEP_DAYS;
        let expected_times = [tt, tt - h, tt + h]
            .map(|at| crate::time::default_timescale_tt_from_ut1(at - 123.0 / 86400.0).unwrap());
        assert_eq!(*times.lock().unwrap(), expected_times);
        times.lock().unwrap().clear();
        let refusing = Rows {
            times: times.clone(),
            refuse: true,
        };
        let error = session
            .calc_with_acquisitions(
                tt,
                MARS,
                flags,
                &PositionAcquisitions {
                    anchor_rows: Some(&refusing),
                    ..Default::default()
                },
            )
            .unwrap_err();
        assert_eq!(error, marker());
        assert_eq!(times.lock().unwrap().len(), 1);
    }
}
