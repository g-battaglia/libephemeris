//! Twenty-component physical phenomena over explicit LEB or local SPK feeds.
use super::*;
use crate::photometry::SealedPhenoServing;
use crate::positions::{MainImageAcquisition, StateProvider};

/// Borrowed public position acquisition for nested phenomena geometry.
/// Arguments are the actually consumed UT1 JD, body and unchanged stage flags.
/// Return all six spherical position/rate components in degrees/AU and per day;
/// only the first three feed photometry. The callback is synchronous and may
/// fail with its original owned Error. It is never retained by a Session or
/// invoked for the direct physical-vector branch.
pub type PhenoPosition<'a> = dyn Fn(f64, i32, i32) -> Result<[f64; 6], Error> + 'a;

/// Same-call owner of an actual nested exception caught by the direct LEB
/// attempt. True consumes that handled failure before original-Time retry;
/// false retains it, and an owned handler error stops before retry.
pub type PhenoPositionRetry<'a> = dyn Fn(&Error) -> Result<bool, Error> + 'a;

/// Original observation Time's consumed UT1 scalar for nested public positions.
/// This is independent of the configured Delta-T model used inside calc_ut.
/// The supplier may retain an opaque observation and propagates its own error.
pub type PhenoObservationUt1<'a> = dyn Fn() -> Result<f64, Error> + 'a;

/// Borrowed native position acquisition, retry classification and UT1 supplier.
/// Positions contain six components in degrees/AU and daily rates. All
/// positional geometry uses the checked native kernel. Nothing is retained
/// beyond this call.
#[derive(Clone, Copy)]
pub struct PhenoCallbacks<'a> {
    position: &'a PhenoPosition<'a>,
    retry: Option<&'a PhenoPositionRetry<'a>>,
    observation_ut1: Option<&'a PhenoObservationUt1<'a>>,
}

impl<'a> PhenoCallbacks<'a> {
    /// Borrow one native UT1 position stage.
    pub fn new(position: &'a PhenoPosition<'a>) -> Self {
        Self {
            position,
            retry: None,
            observation_ut1: None,
        }
    }
    /// Attach the owner of an actually caught direct-attempt exception.
    pub fn with_retry(mut self, retry: &'a PhenoPositionRetry<'a>) -> Self {
        self.retry = Some(retry);
        self
    }
    /// Borrow the actual observation's UT1 descriptor at a point/asteroid
    /// branch. Physical-disc vectors never consume this supplier. Without it,
    /// ordinary checked Core callers retain their existing independent clock.
    pub fn with_observation_ut1(mut self, observation: &'a PhenoObservationUt1<'a>) -> Self {
        self.observation_ut1 = Some(observation);
        self
    }

    /// Borrow the original observation's optional UT1 supplier without reading
    /// its Time descriptor or evaluating any Delta-T/model primitive.
    pub fn observation_ut1(&self) -> Option<&'a PhenoObservationUt1<'a>> {
        self.observation_ut1
    }
    /// Borrow the actual position acquisition, without performing a read.
    pub fn position(&self) -> &'a PhenoPosition<'a> {
        self.position
    }
}

/// Physical phenomena: phase angle (deg), illuminated fraction, elongation
/// (deg), apparent diameter (deg), magnitude and Moon parallax (deg), followed
/// by the remaining producer slots. Preserve all twenty components.
#[derive(Debug, Clone, PartialEq)]
pub struct Phenomena {
    /// Full native-layout attributes, including unused zero slots.
    pub values: [f64; 20],
    /// Normalized compatibility request flags; PHEN acquisition retains its
    /// calculation masks. Physical engine/data identity is separate.
    pub returned_flags: i32,
    /// Session configuration advisories, or empty.
    pub diagnostic: String,
}

impl Phenomena {
    /// Sun-target-observer angle in degrees.
    pub fn phase_angle(&self) -> f64 {
        self.values[0]
    }
    /// Illuminated fraction of the apparent disc.
    pub fn phase(&self) -> f64 {
        self.values[1]
    }
    /// Angular separation from the Sun in degrees.
    pub fn elongation(&self) -> f64 {
        self.values[2]
    }
    /// Apparent diameter in degrees.
    pub fn diameter(&self) -> f64 {
        self.values[3]
    }
    /// Apparent visual magnitude.
    pub fn magnitude(&self) -> f64 {
        self.values[4]
    }
}
impl Session {
    /// Signed difference of target and Sun position component zero at scalar
    /// UT1, with the corresponding evening-star predicate.
    ///
    /// The target is calculated first, then the Sun, using the original epoch
    /// and flags for both calls. The difference is reduced by one turn when
    /// above 180 or below -180. Default spherical output gives degrees;
    /// RADIANS and XYZ retain literal component subtraction and the same
    /// thresholds. A strictly positive result denotes an evening star.
    /// Position acquisition owns timescales, centers, sidereal frames and
    /// observer requirements. Its first error propagates unchanged. The Sun
    /// itself still requires both position calls in this entry.
    pub fn get_elongation_from_sun(
        &self,
        jd_ut: f64,
        body: i32,
        flags: i32,
    ) -> Result<(f64, bool), Error> {
        let target = self.calc_ut(jd_ut, body, flags)?;
        let sun = self.calc_ut(jd_ut, SUN, flags)?;
        let signed = crate::photometry::signed_elongation_deg(target.values[0], sun.values[0]);
        Ok((signed, signed > 0.0))
    }

    /// Signed target-minus-Sun component-zero difference at scalar UT1.
    /// Units, one-turn reduction, flag semantics and acquisition errors follow
    /// [`Self::get_elongation_from_sun`].
    pub fn get_signed_elongation(&self, jd_ut: f64, body: i32, flags: i32) -> Result<f64, Error> {
        self.get_elongation_from_sun(jd_ut, body, flags)
            .map(|(signed, _)| signed)
    }

    /// Whether a target is a morning star at scalar UT1.
    /// The Sun returns false before validating the epoch, flags, observer or
    /// provider. Other targets negate the strictly-positive evening predicate;
    /// zero and unordered differences therefore report true. Acquisition and
    /// component semantics follow [`Self::get_elongation_from_sun`].
    pub fn is_morning_star(&self, jd_ut: f64, body: i32, flags: i32) -> Result<bool, Error> {
        if body == SUN {
            return Ok(false);
        }
        self.get_elongation_from_sun(jd_ut, body, flags)
            .map(|(_, evening)| !evening)
    }

    /// Whether a target is an evening star at scalar UT1.
    /// The Sun returns false without provider, epoch, flag or observer checks.
    /// Other targets require a strictly positive difference, with acquisition
    /// and component semantics from [`Self::get_elongation_from_sun`].
    pub fn is_evening_star(&self, jd_ut: f64, body: i32, flags: i32) -> Result<bool, Error> {
        if body == SUN {
            return Ok(false);
        }
        self.get_elongation_from_sun(jd_ut, body, flags)
            .map(|(_, evening)| evening)
    }

    /// Elongation direction at scalar UT1: "eastern" for a strictly positive
    /// difference and "western" otherwise. The Sun returns "none" without
    /// provider, epoch, flag or observer checks. Other targets retain the
    /// acquisition and component semantics of [`Self::get_elongation_from_sun`].
    pub fn get_elongation_type(
        &self,
        jd_ut: f64,
        body: i32,
        flags: i32,
    ) -> Result<&'static str, Error> {
        if body == SUN {
            return Ok("none");
        }
        self.get_elongation_from_sun(jd_ut, body, flags)
            .map(|(_, evening)| if evening { "eastern" } else { "western" })
    }

    /// Physical phenomena at scalar UT1. The direct sealed route uses the
    /// configured Delta-T model; a stored-reader miss retries the complete
    /// composition with the independent global timescale. Explicit
    /// NOABERR/NOGDEFL select that vector route directly. Range errors retain
    /// the actual support body, queried TT and that body's stored bounds.
    /// NativeJpl instead observes physical-body vectors at the original
    /// global Time; nested public model/asteroid calls keep configured clocks.
    /// Raw split-directory and public callback range contracts remain separate.
    /// All twenty values, including degenerate-body sentinels, are preserved.
    /// Nondegenerate topocentric composition serves through the context
    /// observer; without one it is the typed missing-observer refusal, never an
    /// observer-free substitute.
    /// NativeJpl Earth/ECL_NUT sentinels precede observer and registry acquisition.
    pub fn pheno_ut(&self, jd_ut: f64, body: i32, flags: i32) -> Result<Phenomena, Error> {
        self.pheno_entry(jd_ut, body, flags, true, None, None)
    }
    /// Physical phenomena at scalar TT. The direct route preserves the
    /// configured TT→UT1→TT round trip. Its whole-composition vector retry
    /// instead consumes the original TT and an independent UT1 partner.
    /// NativeJpl directly retains original TT and its independent UT1 partner,
    /// including both parts through target retardation and Topos rotation.
    /// Observer requirements, errors and twenty-slot shape follow `pheno_ut`.
    pub fn pheno(&self, jd_tt: f64, body: i32, flags: i32) -> Result<Phenomena, Error> {
        self.pheno_entry(jd_tt, body, flags, false, None, None)
    }

    /// UT1 phenomena with a call-borrowed public calc_ut stage. Existing time,
    /// observer and provider acquisition precede the callback exactly as in
    /// pheno_ut. Positional bodies consume target then Sun with TRUEPOS/TOPOCTR;
    /// asteroid tables consume target, Sun, then heliocentric target with their
    /// existing masks. Physical-disc vectors do not consume the callback.
    /// All callback errors propagate unchanged; no fallback or resource copy is
    /// added. The caller keeps its callback/resources alive until return.
    pub fn pheno_ut_with_calc_ut(
        &self,
        jd_ut: f64,
        body: i32,
        flags: i32,
        calc_ut: &PhenoPosition<'_>,
    ) -> Result<Phenomena, Error> {
        self.pheno_ut_with_callbacks(jd_ut, body, flags, PhenoCallbacks::new(calc_ut))
    }

    /// TT counterpart of pheno_ut_with_calc_ut. Nested public positions still
    /// consume the original route's actual derived UT1, not the caller TT or a
    /// synthetic epoch. Independent frame/vector clocks and all twenty result
    /// components retain pheno's ordinary semantics.
    pub fn pheno_with_calc_ut(
        &self,
        jd_tt: f64,
        body: i32,
        flags: i32,
        calc_ut: &PhenoPosition<'_>,
    ) -> Result<Phenomena, Error> {
        self.pheno_with_callbacks(jd_tt, body, flags, PhenoCallbacks::new(calc_ut))
    }

    /// Physical phenomena at UT1 with a reusable selected-main image owner.
    /// The owner is consumed only for nondegenerate NativeJpl calculations;
    /// finite-epoch validation and Earth/ECL_NUT sentinels precede acquisition.
    /// Every geometry sample shares the acquired immutable image and current
    /// observer/time settings. Other provider modes retain their ordinary route.
    pub fn pheno_ut_with_main_image(
        &self,
        jd_ut: f64,
        body: i32,
        flags: i32,
        main_image: &dyn MainImageAcquisition,
    ) -> Result<Phenomena, Error> {
        self.pheno_entry(jd_ut, body, flags, true, None, Some(main_image))
    }

    /// TT counterpart of [`Self::pheno_ut_with_main_image`], preserving its
    /// original TT and independent UT1 partner throughout the photometry.
    /// Main-image acquisition errors, observer errors and physical range
    /// failures propagate without an observer-free or coefficient fallback.
    pub fn pheno_with_main_image(
        &self,
        jd_tt: f64,
        body: i32,
        flags: i32,
        main_image: &dyn MainImageAcquisition,
    ) -> Result<Phenomena, Error> {
        self.pheno_entry(jd_tt, body, flags, false, None, Some(main_image))
    }

    /// UT1 public-position counterpart with a fallible direct-LEB catch owner.
    /// Only a caught callback failure consumes retry; its actual class and
    /// exception lifetime belong to the supplied handler. Successful recovery
    /// uses original UT1 and the existing independent vector clocks. Native
    /// and explicit vector requests never invoke the handler. Physical-vector
    /// formulas, typed reader catches and all twenty slots remain unchanged.
    pub fn pheno_ut_with_calc_ut_retry(
        &self,
        jd_ut: f64,
        body: i32,
        flags: i32,
        calc_ut: &PhenoPosition<'_>,
        retry: &PhenoPositionRetry<'_>,
    ) -> Result<Phenomena, Error> {
        self.pheno_ut_with_callbacks(
            jd_ut,
            body,
            flags,
            PhenoCallbacks::new(calc_ut).with_retry(retry),
        )
    }

    /// TT public-position counterpart with the same caught-exception owner.
    /// The first direct stage consumes configured-model UT1; recovery consumes
    /// the independent UT1 of the original caller TT. Neither clock is replaced
    /// by the first attempted stage's epoch.
    pub fn pheno_with_calc_ut_retry(
        &self,
        jd_tt: f64,
        body: i32,
        flags: i32,
        calc_ut: &PhenoPosition<'_>,
        retry: &PhenoPositionRetry<'_>,
    ) -> Result<Phenomena, Error> {
        self.pheno_with_callbacks(
            jd_tt,
            body,
            flags,
            PhenoCallbacks::new(calc_ut).with_retry(retry),
        )
    }

    /// UT1 phenomena with borrowed native positions and retry classification.
    /// All public acquisitions precede geometry. Actual arithmetic refusals
    /// participate in the same direct-attempt catch/retry lifecycle; successful
    /// derived NaNs keep the source comparison branches instead of a finite
    /// state constructor. The ordinary physical-vector route is unchanged.
    pub fn pheno_ut_with_callbacks(
        &self,
        jd_ut: f64,
        body: i32,
        flags: i32,
        callbacks: PhenoCallbacks<'_>,
    ) -> Result<Phenomena, Error> {
        self.pheno_entry(jd_ut, body, flags, true, Some(callbacks), None)
    }

    /// TT counterpart with original configured/independent retry clocks.
    /// Scalar arithmetic belongs only to the completed public positional
    /// geometry, never the physical vector or original-Time acquisition.
    pub fn pheno_with_callbacks(
        &self,
        jd_tt: f64,
        body: i32,
        flags: i32,
        callbacks: PhenoCallbacks<'_>,
    ) -> Result<Phenomena, Error> {
        self.pheno_entry(jd_tt, body, flags, false, Some(callbacks), None)
    }
    /// Bind one temporary photometry serving view; no result/image cache copies.
    #[allow(clippy::too_many_arguments)]
    fn pheno_entry(
        &self,
        jd: f64,
        body: i32,
        flags: i32,
        is_ut: bool,
        callbacks: Option<PhenoCallbacks<'_>>,
        main_image: Option<&dyn MainImageAcquisition>,
    ) -> Result<Phenomena, Error> {
        if !jd.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd });
        }
        // Phenomena owns its acquisition flags. In particular, the position
        // normalizer implies NOABERR/NOGDEFL for TRUEPOS; those implied bits
        // must not become an explicit request for the non-LEB phenomena leg.
        let calculation_flags = flags;
        let flags = crate::positions::normalize_calc_flags(flags, false);
        // Topocentric composition serves through the context observer
        // without one it keeps the typed
        // missing-observer refusal from the serving layer, never a
        // silent geocentric substitute.
        let site = self
            .ctx
            .observer()
            .map(|o| (o.lon_deg(), o.lat_deg(), o.height_m()));
        if self.ctx.mode() == ProviderMode::NativeJpl {
            // Degenerate Earth/nutation tuples consume no observer or ephemeris.
            // Explicitly supplied complete registries also remain independent
            // of path acquisition; partial target maps still need Earth/Sun.
            let degenerate = matches!(body, EARTH | ECL_NUT);
            let needs_main = if degenerate {
                false
            } else if main_image.is_some() {
                true
            } else {
                let provider = crate::positions::LocalSpkStateProvider::new(
                    self.ctx.spk_registry(),
                    self.ctx.mode(),
                )?;
                let physical = crate::positions::remap_ast_offset(body);
                !provider.has_body(EARTH)
                    || !provider.has_body(SUN)
                    || ((SUN..=PLUTO).contains(&physical) && !provider.has_body(physical))
            };
            let acquired = if needs_main {
                let owner = match main_image {
                    Some(main_image) => main_image.acquire(&self.ctx)?,
                    None => MainImageLifecycle::new().acquire(&self.ctx)?,
                };
                Some(owner.with_current_context(self.ctx.clone())?)
            } else {
                None
            };
            let selected = acquired
                .as_ref()
                .map(|owner| owner.native_provider_context())
                .transpose()?;
            let native = selected.as_ref().unwrap_or(&self.ctx);
            let values = crate::positions::native_jpl_pheno_with_callbacks(
                native,
                self.delta_t_model(),
                self.ctx.effective_iers_delta_t(),
                jd,
                is_ut,
                body,
                calculation_flags,
                callbacks.as_ref(),
            )?;
            return Ok(Phenomena {
                values,
                returned_flags: flags,
                diagnostic: self.diagnostics.join("; "),
            });
        }
        let values = self.with_serving_stack(|stack| {
            let delta = |ut| self.deltat(ut);
            let position = |ut, body, calc_flags| match callbacks {
                Some(callbacks) => (callbacks.position)(ut, body, calc_flags),
                None => self.calc_ut(ut, body, calc_flags).map(|p| p.values),
            };
            let mut serving = SealedPhenoServing::new(stack.serving, &delta, calculation_flags)
                .with_calc_ut(&position);
            if let Some(callbacks) = callbacks {
                if let Some(retry) = callbacks.retry {
                    serving = serving.with_position_retry(retry);
                }
            }
            if let Some(site) = site {
                serving = serving.with_site(site);
            }
            if is_ut {
                serving.pheno_ut(body, jd)
            } else {
                serving.pheno(body, jd)
            }
        })?;
        Ok(Phenomena {
            values,
            returned_flags: flags,
            diagnostic: self.diagnostics.join("; "),
        })
    }
}

#[cfg(test)]
mod elongation_tests {
    use super::*;

    #[test]
    fn sun_predicates_precede_epoch_flags_observer_and_provider() {
        for mode in [
            ProviderMode::SealedLeb,
            ProviderMode::Auto,
            ProviderMode::NativeJpl,
            ProviderMode::HorizonsRemote,
        ] {
            let session = Session::from_context(EphemerisContext::new().with_mode(mode));
            let flags = FLG_TOPOCTR | FLG_BARYCTR | FLG_MOSEPH;
            for epoch in [2451545.125, f64::NAN, f64::INFINITY] {
                assert!(!session.is_morning_star(epoch, SUN, flags).unwrap());
                assert!(!session.is_evening_star(epoch, SUN, flags).unwrap());
                assert_eq!(
                    session.get_elongation_type(epoch, SUN, flags).unwrap(),
                    "none"
                );
            }
        }
    }

    #[test]
    fn signed_sun_and_non_sun_entries_retain_position_errors() {
        let session =
            Session::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
        assert!(session.get_signed_elongation(2451545.125, SUN, 0).is_err());
        assert!(session.is_morning_star(2451545.125, MARS, 0).is_err());
        for body in [SUN, MARS] {
            let error = session
                .get_signed_elongation(f64::INFINITY, body, 0)
                .unwrap_err();
            assert!(matches!(
                error.primary_error(),
                Error::NonFiniteArgument { value } if *value == f64::INFINITY
            ));
        }
    }
}

#[cfg(test)]
mod callback_tests {
    use super::*;
    use std::cell::RefCell;

    /// Empty index metadata admits a serving view but has no astronomical
    /// coefficients. These controls prove dispatch/ownership, not science.
    fn callback_session() -> Session {
        callback_session_with_aliases(&[])
    }

    /// Optional zero-valued toy aliases satisfy only the unchanged vector
    /// admission. The asteroid callback branch consumes no coefficient words.
    fn callback_session_with_aliases(aliases: &[i32]) -> Session {
        use crate::leb::{
            COORD_ICRS_BARY, HEADER_SIZE, LEB1_BODY_ENTRY_SIZE, LebFile, MAGIC_LEB1,
            SECTION_DIR_SIZE,
        };
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&MAGIC_LEB1);
        for value in [1u32, 1, aliases.len() as u32] {
            bytes.extend_from_slice(&value.to_le_bytes());
        }
        for value in [2450000.0f64, 2470000.0, 0.0] {
            bytes.extend_from_slice(&value.to_le_bytes());
        }
        bytes.extend_from_slice(&0u32.to_le_bytes());
        bytes.extend_from_slice(&[0u8; 20]);
        assert_eq!(bytes.len(), HEADER_SIZE);
        bytes.extend_from_slice(&0u32.to_le_bytes());
        bytes.extend_from_slice(&0u32.to_le_bytes());
        bytes.extend_from_slice(&((HEADER_SIZE + SECTION_DIR_SIZE) as u64).to_le_bytes());
        bytes.extend_from_slice(&((aliases.len() * LEB1_BODY_ENTRY_SIZE) as u64).to_le_bytes());
        let data_at = HEADER_SIZE + SECTION_DIR_SIZE + aliases.len() * LEB1_BODY_ENTRY_SIZE;
        for (index, &body) in aliases.iter().enumerate() {
            bytes.extend_from_slice(&body.to_le_bytes());
            bytes.extend_from_slice(&COORD_ICRS_BARY.to_le_bytes());
            bytes.extend_from_slice(&1u32.to_le_bytes());
            for value in [2450000.0f64, 2470000.0, 20000.0] {
                bytes.extend_from_slice(&value.to_le_bytes());
            }
            bytes.extend_from_slice(&0u32.to_le_bytes()); // degree zero
            bytes.extend_from_slice(&3u32.to_le_bytes()); // three components
            bytes.extend_from_slice(&((data_at + index * 24) as u64).to_le_bytes());
        }
        bytes.resize(data_at + aliases.len() * 24, 0);
        let mut context = EphemerisContext::new()
            .with_mode(ProviderMode::SealedLeb)
            .with_user_delta_t(Some(0.25))
            .with_observer(Observer::new(13.0, 52.0, 34.0).unwrap());
        context.attach_leb_file(
            LebFile::from_bytes("mem://pheno-callback-empty-index.leb1".into(), bytes).unwrap(),
        );
        Session::from_context(context)
    }

    fn position(body: i32) -> [f64; 6] {
        // Unconsumed rates must not be preflighted or used in photometry.
        [
            if body == SUN { 0.0 } else { 90.0 },
            0.0,
            1.0,
            f64::NAN,
            f64::INFINITY,
            f64::NEG_INFINITY,
        ]
    }

    #[test]
    fn pheno_callback_keeps_ut1_two_stage_mask_and_all_output_slots() {
        let session = callback_session();
        let flags = FLG_SWIEPH
            | FLG_TRUEPOS
            | FLG_TOPOCTR
            | FLG_NONUT
            | FLG_HELCTR
            | FLG_SPEED
            | FLG_XYZ
            | FLG_RADIANS
            | FLG_EQUATORIAL
            | FLG_SIDEREAL;
        let calls = RefCell::new(Vec::new());
        let callback = |ut, body, stage_flags| {
            calls.borrow_mut().push((ut, body, stage_flags));
            Ok(position(body))
        };
        let result = session
            .pheno_ut_with_calc_ut(2460400.25, CUPIDO, flags, &callback)
            .unwrap();
        assert_eq!(
            *calls.borrow(),
            [
                (2460400.25, CUPIDO, FLG_TRUEPOS | FLG_TOPOCTR),
                (2460400.25, SUN, FLG_TRUEPOS | FLG_TOPOCTR),
            ]
        );
        assert!((result.phase_angle() - 45.0).abs() < 1e-12);
        assert!((result.elongation() - 90.0).abs() < 1e-12);
        assert!(result.values[3..].iter().all(|word| word.to_bits() == 0));
        assert_eq!(
            result.returned_flags,
            crate::positions::normalize_calc_flags(flags, false)
        );
    }

    #[test]
    fn pheno_callback_tt_passes_actual_configured_ut1() {
        let session = callback_session();
        let calls = RefCell::new(Vec::new());
        let callback = |ut, body, stage_flags| {
            calls.borrow_mut().push((ut, body, stage_flags));
            Ok(position(body))
        };
        session
            .pheno_with_calc_ut(2460400.25, CUPIDO, FLG_SWIEPH | FLG_NONUT, &callback)
            .unwrap();
        assert_eq!(
            *calls.borrow(),
            [(2460400.0, CUPIDO, 0), (2460400.0, SUN, 0)]
        );
    }

    #[test]
    fn pheno_callback_asteroid_stops_at_each_original_fallible_stage() {
        let session = callback_session_with_aliases(&[EARTH, SUN]);
        // Source NOABERR bypasses the direct LEB attempt. The ensuing vector
        // route admits Earth/Sun before the public asteroid position stages.
        let flags = FLG_NOABERR | FLG_TRUEPOS | FLG_NONUT | FLG_TOPOCTR | FLG_SPEED | FLG_XYZ;
        let expected = [
            (CERES, FLG_TRUEPOS | FLG_NONUT | FLG_TOPOCTR),
            (SUN, FLG_TRUEPOS | FLG_NONUT | FLG_TOPOCTR),
            (CERES, FLG_TRUEPOS | FLG_NONUT | FLG_HELCTR),
        ];
        for refused in 0..3 {
            let calls = RefCell::new(Vec::new());
            let callback = |ut: f64, body, stage_flags| {
                assert_eq!(ut.to_bits(), 2460400.25f64.to_bits());
                let mut seen = calls.borrow_mut();
                let ordinal = seen.len();
                seen.push((body, stage_flags));
                if ordinal == refused {
                    Err(Error::ConfigurationMissing {
                        missing_config: "actual_callback_owner",
                        suggestion: "actual_callback_refusal",
                    })
                } else {
                    Ok(position(body))
                }
            };
            let error = session
                .pheno_ut_with_calc_ut(2460400.25, CERES, flags, &callback)
                .unwrap_err();
            assert!(matches!(
                error,
                Error::ConfigurationMissing {
                    missing_config: "actual_callback_owner",
                    suggestion: "actual_callback_refusal"
                }
            ));
            assert_eq!(*calls.borrow(), expected[..=refused]);
        }
    }

    #[test]
    fn pheno_callback_preserves_earlier_guards_and_physical_vector_route() {
        let never = |_: f64, _: i32, _: i32| -> Result<[f64; 6], Error> {
            panic!("unused positional callback")
        };
        let session = callback_session();
        assert!(matches!(
            session.pheno_ut_with_calc_ut(f64::INFINITY, CUPIDO, 0, &never),
            Err(Error::NonFiniteArgument { value }) if value == f64::INFINITY
        ));
        let empty = Session::from_context(EphemerisContext::new());
        assert!(
            empty
                .pheno_ut_with_calc_ut(2460400.25, CUPIDO, 0, &never)
                .is_err()
        );
        // The empty index forces the unchanged physical provider failure.
        // It must never turn this direct Sun request into public positions.
        let ordinary = session.pheno_ut(2460400.25, SUN, FLG_TRUEPOS);
        let callback = session.pheno_ut_with_calc_ut(2460400.25, SUN, FLG_TRUEPOS, &never);
        assert_eq!(format!("{ordinary:?}"), format!("{callback:?}"));
        assert!(ordinary.is_err());
    }

    #[test]
    fn pheno_callback_retry_keeps_original_ut1_and_distinct_tt_clocks() {
        let session = callback_session_with_aliases(&[EARTH, SUN]);
        let jd = 2460400.25;
        for is_ut in [false, true] {
            let calls = RefCell::new(Vec::new());
            let retried = RefCell::new(0usize);
            let callback = |ut, body, flags| {
                let mut calls = calls.borrow_mut();
                calls.push((ut, body, flags));
                if calls.len() == 1 {
                    Err(Error::ConfigurationMissing {
                        missing_config: "caught_callback",
                        suggestion: "handle the first actual exception",
                    })
                } else {
                    Ok(position(body))
                }
            };
            let retry = |error: &Error| {
                assert!(matches!(
                    error,
                    Error::ConfigurationMissing {
                        missing_config: "caught_callback",
                        ..
                    }
                ));
                *retried.borrow_mut() += 1;
                Ok(true)
            };
            let result = if is_ut {
                session.pheno_ut_with_calc_ut_retry(jd, CUPIDO, FLG_SWIEPH, &callback, &retry)
            } else {
                session.pheno_with_calc_ut_retry(jd, CUPIDO, FLG_SWIEPH, &callback, &retry)
            }
            .unwrap();
            let first = if is_ut { jd } else { jd - 0.25 };
            let second = if is_ut {
                jd
            } else {
                crate::time::default_timescale_ut1_from_tt(jd).unwrap()
            };
            assert_eq!(
                *calls.borrow(),
                [(first, CUPIDO, 0), (second, CUPIDO, 0), (second, SUN, 0)]
            );
            assert_eq!(*retried.borrow(), 1);
            if !is_ut {
                assert_ne!(first.to_bits(), second.to_bits());
            }
            assert!((result.phase_angle() - 45.0).abs() < 1e-12);
        }
    }

    #[test]
    fn pheno_callback_final_retry_error_is_not_caught_again() {
        let session = callback_session_with_aliases(&[EARTH, SUN]);
        let calls = RefCell::new(0usize);
        let handled = RefCell::new(0usize);
        let callback = |_: f64, _: i32, _: i32| {
            *calls.borrow_mut() += 1;
            Err(Error::ConfigurationMissing {
                missing_config: if *calls.borrow() == 1 {
                    "first"
                } else {
                    "final"
                },
                suggestion: "actual owning callback",
            })
        };
        let retry = |_: &Error| {
            *handled.borrow_mut() += 1;
            Ok(true)
        };
        assert!(matches!(
            session.pheno_ut_with_calc_ut_retry(2460400.25, CUPIDO, 0, &callback, &retry),
            Err(Error::ConfigurationMissing {
                missing_config: "final",
                ..
            })
        ));
        assert_eq!(*calls.borrow(), 2);
        assert_eq!(*handled.borrow(), 1);
    }

    #[test]
    fn pheno_callback_retry_handler_failure_precedes_vector_acquisition() {
        let session = callback_session_with_aliases(&[EARTH, SUN]);
        let calls = RefCell::new(0usize);
        let callback = |_: f64, _: i32, _: i32| {
            *calls.borrow_mut() += 1;
            Err(Error::ConfigurationMissing {
                missing_config: "first",
                suggestion: "actual owning callback",
            })
        };
        let retry = |_: &Error| {
            Err(Error::ConfigurationMissing {
                missing_config: "handler",
                suggestion: "actual owning retry handler",
            })
        };
        assert!(matches!(
            session.pheno_ut_with_calc_ut_retry(2460400.25, CUPIDO, 0, &callback, &retry),
            Err(Error::ConfigurationMissing {
                missing_config: "handler",
                ..
            })
        ));
        assert_eq!(*calls.borrow(), 1);
    }

    #[test]
    fn pheno_callback_refuses_nonfinite_positions_after_acquisition() {
        let session = callback_session();
        for is_ut in [false, true] {
            let calls = RefCell::new(Vec::new());
            let callback = |ut, body, flags| {
                calls.borrow_mut().push((ut, body, flags));
                Ok(if body == SUN {
                    [0.0, 0.0, 1.0, f64::INFINITY, f64::NAN, f64::NEG_INFINITY]
                } else {
                    [
                        f64::NAN,
                        0.0,
                        1.0,
                        f64::NAN,
                        f64::INFINITY,
                        f64::NEG_INFINITY,
                    ]
                })
            };
            let error = if is_ut {
                session.pheno_ut_with_calc_ut(2460400.25, CUPIDO, FLG_SWIEPH, &callback)
            } else {
                session.pheno_with_calc_ut(2460400.25, CUPIDO, FLG_SWIEPH, &callback)
            }
            .unwrap_err();
            let epoch = if is_ut { 2460400.25 } else { 2460400.0 };
            assert_eq!(*calls.borrow(), [(epoch, CUPIDO, 0), (epoch, SUN, 0)]);
            assert!(matches!(error, Error::NonFiniteArgument { value } if value.is_nan()));
        }
    }
}
