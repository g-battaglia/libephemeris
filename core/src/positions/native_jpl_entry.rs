//! Native JPL public-entry acquisition and call-scoped callback assembly.
//!
//! Refreshed minor-body registries and main-image owners remain alive through
//! the calculation callback. Entry metadata preserves original clocks, fallback
//! flags and trace scopes; shared geometry and range policies stay in routing.

#[path = "native_jpl_entry/hypothetical_models.rs"]
mod hypothetical_models;
#[path = "native_jpl_entry/minor_fallback.rs"]
mod minor_fallback;

use super::{
    LocalSpkStateProvider, NativeTraceInvocation, PublicModelDirectory, local_jpl_lunar_position,
    local_jpl_minor_body_candidate, local_jpl_minor_body_servable, local_jpl_position,
    local_jpl_position_with_public_entry, local_jpl_sidereal_word, local_jpl_target,
    minor_requires_main_image, model_serving, native_sidereal_frame_request, native_source_result,
    publish_native_trace_source, validate_local_jpl_public_epoch, vector_self_observation,
    with_native_trace_body,
};
use crate::constants::{
    CUPIDO, EARTH, FLG_BARYCTR, FLG_HELCTR, INTP_APOG, INTP_PERG, MEAN_APOG, MEAN_NODE, OSCU_APOG,
    PLUTO, SUN, TRUE_NODE, WALDEMATH,
};
use crate::context::ProviderMode;
use crate::errors::{Error, Result};
use crate::positions::{
    AnalyticalModelCallback, CalcEngine, CartesianCompute, JplLeg, PositionResult, ServingRequest,
    calc_ut_with_cartesian, calc_with_cartesian, consume_non_echoed_flags, ecl_nut_result,
    ecl_nut_result_ut, forbid_jpl_in_sealed_mode, normalize_calc_flags,
    validate_barycentric_moseph,
};
use crate::time::deltat_with_config;

/// Calculate with an optional native analytical-model provider.
///
/// The epoch, output units, flags and errors follow
/// [`crate::positions::calc_native_jpl_entry`]. The borrowed model callback is
/// consumed synchronously only when the selected body needs an analytical state.
#[allow(clippy::too_many_arguments)]
pub fn calc_native_jpl_entry_with_models(
    models: Option<&AnalyticalModelCallback<'_>>,
    req: ServingRequest<'_>,
    epoch: f64,
    body: i32,
    flags: i32,
    ut_entry: bool,
) -> Result<(PositionResult, i32)> {
    calc_native_jpl_entry_with_main_image(
        models, None, None, req, epoch, body, flags, ut_entry, None,
    )
}

/// Calculate with optional native analytical models and dynamic sidereal anchors.
///
/// The epoch, output units, flags and errors follow
/// [`crate::positions::calc_native_jpl_entry`]. Providers are borrowed
/// synchronously; the calculation does not retain their callbacks.
#[allow(clippy::too_many_arguments)]
pub fn calc_native_jpl_entry_with_acquisitions(
    models: Option<&AnalyticalModelCallback<'_>>,
    live_anchors: Option<&dyn crate::sidereal::DynamicAnchorProvider>,
    req: ServingRequest<'_>,
    epoch: f64,
    body: i32,
    flags: i32,
    ut_entry: bool,
) -> Result<(PositionResult, i32)> {
    calc_native_jpl_entry_with_main_image(
        models,
        live_anchors,
        None,
        req,
        epoch,
        body,
        flags,
        ut_entry,
        None,
    )
}

/// Calculate with the registered-SPK, auto-SPK, ASSIST and analytical minor-body chain.
///
/// The epoch, output units and flag handling follow
/// [`crate::positions::calc_native_jpl_entry`]. The context selects the precision
/// and network policy; missing data and coverage failures retain that policy's
/// errors. The numerical supplier reuses its owned ephemerides across calls.
#[allow(clippy::too_many_arguments)]
pub fn calc_native_jpl_entry_with_minor_assist(
    models: Option<&AnalyticalModelCallback<'_>>,
    live_anchors: Option<&dyn crate::sidereal::DynamicAnchorProvider>,
    minor_assist: Option<&dyn crate::engine_api::MinorAssistAcquisition>,
    req: ServingRequest<'_>,
    epoch: f64,
    body: i32,
    flags: i32,
    ut_entry: bool,
) -> Result<(PositionResult, i32)> {
    calc_native_jpl_entry_with_main_image(
        models,
        live_anchors,
        minor_assist,
        req,
        epoch,
        body,
        flags,
        ut_entry,
        None,
    )
}

/// Calculate with a reusable main-image owner acquired only when geometry needs it.
///
/// `epoch` is a UT1 Julian day when `ut_entry` is true and a TT Julian day
/// otherwise. Coordinates, rates and returned flags follow
/// [`crate::positions::calc_native_jpl_entry`]. Actual body/model coverage still
/// bounds the calculation; retaining an image does not extend its domain.
///
/// The configuration snapshot remains immutable. Refreshed target registrations
/// and the acquired main image stay owned until calculation and tracing finish.
/// Invalid flags, nonfinite time, data/coverage failures and numerical errors
/// propagate through the selected precision and network policy.
#[allow(clippy::too_many_arguments)]
pub fn calc_native_jpl_entry_with_main_image(
    models: Option<&AnalyticalModelCallback<'_>>,
    live_anchors: Option<&dyn crate::sidereal::DynamicAnchorProvider>,
    minor_assist: Option<&dyn crate::engine_api::MinorAssistAcquisition>,
    req: ServingRequest<'_>,
    epoch: f64,
    body: i32,
    flags: i32,
    ut_entry: bool,
    main_image: Option<&dyn super::MainImageAcquisition>,
) -> Result<(PositionResult, i32)> {
    calc_native_jpl_entry_with_owned_minor(
        models,
        live_anchors,
        minor_assist,
        None,
        req,
        epoch,
        body,
        flags,
        ut_entry,
        main_image,
    )
}

impl crate::engine_api::NativeMinorAssist<'_> {
    /// Concrete default-entry ownership, never inferred from a public callback.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn calc_native_position_with_main_image(
        &self,
        models: Option<&AnalyticalModelCallback<'_>>,
        live_anchors: Option<&dyn crate::sidereal::DynamicAnchorProvider>,
        req: ServingRequest<'_>,
        epoch: f64,
        body: i32,
        flags: i32,
        ut_entry: bool,
        main_image: Option<&dyn super::MainImageAcquisition>,
    ) -> Result<(PositionResult, i32)> {
        calc_native_jpl_entry_with_owned_minor(
            models,
            live_anchors,
            Some(self),
            Some(self),
            req,
            epoch,
            body,
            flags,
            ut_entry,
            main_image,
        )
    }

    /// Differentiate the actual reception Earth position on its owned records.
    /// Stored Type3 velocities are independent channels and are never used as
    /// position tangents. Samples retain the complete split Circular179 map.
    pub(crate) fn earth_position_tt_rate(
        context: &crate::context::EphemerisContext,
        tt: f64,
        central_position: [f64; 3],
    ) -> Result<[f64; 3]> {
        let registry = context.spk_registry();
        let target = local_jpl_target(&JplLeg::new(registry), EARTH)?;
        let clock = |epoch: f64| -> Result<crate::types::TwoPartEpoch> {
            let whole = epoch.floor();
            let fraction = epoch - whole;
            Ok(crate::types::TwoPartEpoch {
                hi: whole,
                lo: fraction
                    + super::reference_frames::tdb_minus_tt_seconds(whole, fraction)? / 86400.0,
            })
        };
        let anchor = clock(tt)?;
        let bounds = registry.geometric_ssb_record_offsets(target, anchor.hi, anchor.lo)?;
        let anchor_value = anchor.hi + anchor.lo;
        let (start, end) = crate::time::tt_interval_from_tdb(
            (anchor_value + bounds.0).next_up(),
            (anchor_value + bounds.1).next_down(),
        )?;
        let provider = LocalSpkStateProvider::new(registry, context.mode())?;
        crate::positions::derivative::cartesian_rates_on_owned_domain(
            &|offset| {
                if offset == 0.0 {
                    return Ok(Some(central_position));
                }
                let Some(epoch) = represented_predicted_sample_epoch(tt, offset) else {
                    return Ok(None);
                };
                let at = clock(epoch)?;
                let Some(relative) = PredictedSpkReads::offsets(at, anchor) else {
                    return Ok(None);
                };
                if relative.0 < bounds.0 || relative.1 > bounds.1 {
                    return Ok(None);
                }
                provider
                    .nodaps_state_two_part(
                        EARTH,
                        crate::types::TwoPartEpoch {
                            hi: epoch.floor(),
                            lo: epoch - epoch.floor(),
                        },
                    )
                    .map(|state| Some(state.pos()))
            },
            1.0 / 64.0,
            (start - tt, end - tt),
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn calc_native_jpl_entry_with_owned_minor(
    models: Option<&AnalyticalModelCallback<'_>>,
    live_anchors: Option<&dyn crate::sidereal::DynamicAnchorProvider>,
    minor_assist: Option<&dyn crate::engine_api::MinorAssistAcquisition>,
    owned_minor: Option<&crate::engine_api::NativeMinorAssist<'_>>,
    req: ServingRequest<'_>,
    epoch: f64,
    body: i32,
    flags: i32,
    ut_entry: bool,
    main_image: Option<&dyn super::MainImageAcquisition>,
) -> Result<(PositionResult, i32)> {
    validate_barycentric_moseph(flags, if ut_entry { "calc_ut" } else { "calc" })?;
    if !epoch.is_finite() {
        return Err(Error::NonFiniteArgument { value: epoch });
    }
    if main_image.is_none()
        && ut_entry
        && matches!(body, SUN..=PLUTO | EARTH)
        && local_jpl_target(&JplLeg::new(req.ctx.spk_registry()), body).is_ok()
    {
        validate_local_jpl_public_epoch(&JplLeg::new(req.ctx.spk_registry()), epoch, body, true)?;
    }
    // This entry is reached only after the stored-data path has declined.
    // Provision a minor target lazily, preserving the original call's settings
    // and retaining only the refreshed immutable SPK resources.
    let refreshed = if local_jpl_minor_body_candidate(body) {
        let tt = if ut_entry {
            crate::smh2016::source_timescale_tt_from_ut1(epoch)?
        } else {
            epoch
        };
        if !local_jpl_minor_body_servable(&JplLeg::new(req.ctx.spk_registry()), body, tt, flags) {
            minor_assist
                .map(|supplier| supplier.refreshed_spk(body, tt, req.ctx))
                .transpose()?
                .flatten()
        } else {
            None
        }
    } else {
        None
    };
    let refreshed_context = refreshed.map(|registry| {
        let mut context = req.ctx.clone();
        context.spk_registry_mut().merge_shared_from(&registry);
        context
    });
    let req = ServingRequest {
        ctx: refreshed_context.as_ref().unwrap_or(req.ctx),
        ..req
    };
    let needs_minor_observer = if minor_assist.is_some() && local_jpl_minor_body_candidate(body) {
        let tt = if ut_entry {
            crate::smh2016::source_timescale_tt_from_ut1(epoch)?
        } else {
            epoch
        };
        minor_requires_main_image(&JplLeg::new(req.ctx.spk_registry()), body, tt, flags)?
    } else {
        false
    };
    let acquired = if req.ctx.mode() == ProviderMode::NativeJpl
        && (matches!(body, SUN..=PLUTO | EARTH | CUPIDO..=WALDEMATH)
            || crate::lunar::is_lunar_point(body)
            || needs_minor_observer)
    {
        if let Some(main_image) = main_image {
            let owner = main_image
                .acquire(req.ctx)?
                .with_current_context(req.ctx.clone())?;
            Some(super::main_image::AcquiredMainImage {
                context: owner.native_provider_context()?,
                fallback_active: owner.fallback_active(),
                main_owner: Some(owner),
            })
        } else if local_jpl_target(&JplLeg::new(req.ctx.spk_registry()), SUN).is_err() {
            Some(super::main_image::acquire_missing_main_image(req.ctx)?)
        } else {
            None
        }
    } else {
        None
    };
    let fallback_active = acquired.as_ref().is_some_and(|image| {
        image
            .main_owner
            .as_ref()
            .map_or(image.fallback_active, |owner| owner.fallback_active())
    });
    let req = ServingRequest {
        ctx: acquired.as_ref().map_or(req.ctx, |image| &image.context),
        ..req
    };
    if ut_entry && matches!(body, SUN..=PLUTO | EARTH) {
        forbid_jpl_in_sealed_mode(req.ctx.mode())?;
        validate_local_jpl_public_epoch(&JplLeg::new(req.ctx.spk_registry()), epoch, body, true)?;
    }
    let _trace_scope = NativeTraceInvocation::enter(Some(body));
    with_local_jpl_entry_outputs_and_fallback_and_models(
        models,
        live_anchors,
        minor_assist,
        req,
        LocalJplEntry {
            original_ut: ut_entry.then_some(epoch),
            fallback_active,
            public_entry: true,
            owned_minor,
        },
        |engine, cartesian| {
            if ut_entry && (CUPIDO..=WALDEMATH).contains(&body) {
                // Retain the original validated UT1 observation while evaluating the derived model epoch.
                let raw = consume_non_echoed_flags(flags, body);
                let norm = normalize_calc_flags(raw, engine.jpl_fallback_active);
                with_native_trace_body(body, body, || {
                    let tt = crate::smh2016::source_timescale_tt_from_ut1(epoch)?;
                    let result =
                        super::serving::calc_tt_body_source_time(engine, tt, body, raw, norm)?;
                    publish_native_trace_source(body)?;
                    Ok(result)
                })
            } else if ut_entry {
                calc_ut_with_cartesian(engine, cartesian, epoch, body, flags)
            } else {
                calc_with_cartesian(engine, cartesian, epoch, body, flags)
            }
        },
    )
}

/// Already acquired scalar entry facts. Grouping these values does not select
/// a callback, provider, frame clock or error constructor.
struct LocalJplEntry<'a> {
    original_ut: Option<f64>,
    fallback_active: bool,
    public_entry: bool,
    owned_minor: Option<&'a crate::engine_api::NativeMinorAssist<'a>>,
}

/// Per-call model-acquisition variant retaining the corresponding static
/// provider/clock/flag pipeline. The borrowed callback is consumed synchronously
/// only at analytical geometric state reads; other bodies never call it.
/// Model callback assembly with an optional genuine scalar range constructor.
/// Numerical/model callbacks and the legacy None route are unchanged.
fn with_local_jpl_entry_outputs_and_fallback_and_models<R>(
    models: Option<&AnalyticalModelCallback<'_>>,
    live_anchors: Option<&dyn crate::sidereal::DynamicAnchorProvider>,
    minor_assist: Option<&dyn crate::engine_api::MinorAssistAcquisition>,
    req: ServingRequest<'_>,
    entry: LocalJplEntry<'_>,
    f: impl FnOnce(&CalcEngine<'_>, &CartesianCompute<'_>) -> Result<R>,
) -> Result<R> {
    let LocalJplEntry {
        original_ut,
        fallback_active,
        public_entry,
        owned_minor,
    } = entry;
    forbid_jpl_in_sealed_mode(req.ctx.mode())?;
    let leg = JplLeg::new(req.ctx.spk_registry());
    let directory = PublicModelDirectory {
        leg: &leg,
        raw: std::cell::RefCell::new(None),
    };
    let site = req
        .observer
        .or_else(|| req.ctx.observer())
        .map(|obs| (obs.lon_deg(), obs.lat_deg(), obs.height_m()));
    let states = LocalSpkStateProvider::new(req.ctx.spk_registry(), req.ctx.mode())?;
    let anchors = crate::stars::SealedAnchorProvider::new(&states);
    let anchors = if let Some(live) = live_anchors {
        anchors.with_dynamic_override(live)
    } else {
        anchors
    };
    let (selection, _) = crate::houses::sidereal_selection_for_context(req.ctx);
    let delta_cfg = req.ctx.delta_t_config(req.delta_t_model, req.iers_enabled);
    let model_deltat = |ut| deltat_with_config(ut, &delta_cfg);
    let ayanamsha =
        crate::sidereal::SealedAyanamshaServing::new(selection, &model_deltat, Some(&anchors));
    let sbdb_cache = std::cell::RefCell::new(crate::minor_bodies::SbdbElementsCache::new());
    let sbdb_elements = |number| {
        crate::minor_bodies::fetch_sbdb_elements(
            number,
            30.0,
            req.ctx.effective_network_policy(),
            &mut sbdb_cache.borrow_mut(),
        )
    };
    let compute_raw = |tt: f64, body: i32, flags: i32| {
        // Abstract lunar points have no heliocentric/barycentric place.
        if matches!(
            body,
            MEAN_NODE | TRUE_NODE | MEAN_APOG | OSCU_APOG | INTP_APOG | INTP_PERG
        ) && flags & (FLG_HELCTR | FLG_BARYCTR) != 0
        {
            if !tt.is_finite() {
                return Err(Error::NonFiniteArgument { value: tt });
            }
            return native_source_result(Ok(([0.0; 6], flags)), "Analytical");
        }
        if crate::lunar::is_lunar_point(body) {
            return native_source_result(
                local_jpl_lunar_position(&leg, tt, original_ut, body, flags, &|ut| {
                    crate::planetary_moons::moon_spk_ayanamsha(ut, flags, &ayanamsha)
                }),
                if matches!(body, TRUE_NODE | OSCU_APOG) {
                    "Skyfield"
                } else {
                    "Analytical"
                },
            );
        }
        // This borrowed model view performs no support reads until its original
        // branch is selected. Publication stays here, outside family reduction.
        if let Some(output) = (hypothetical_models::HypotheticalModels {
            models,
            directory: &directory,
            original_ut,
            model_deltat: &model_deltat,
            ayanamsha: &ayanamsha,
            selection,
        })
        .position(tt, body, flags, site)?
        {
            return native_source_result(Ok((output, flags)), "Analytical");
        }
        if local_jpl_minor_body_candidate(body)
            && !local_jpl_minor_body_servable(&leg, body, tt, flags)
        {
            // Registered/refreshed SPK and automatic provisioning have already
            // had priority. Borrow the remaining policy inputs only on this miss;
            // the entry still owns acquisition scopes and the shared SBDB cache.
            return minor_fallback::MinorBodyFallback {
                leg: &leg,
                context: req.ctx,
                minor_assist,
                owned_minor,
                original_ut,
                ayanamsha_at: &|ut| {
                    crate::planetary_moons::moon_spk_ayanamsha(ut, flags, &ayanamsha)
                },
                sbdb_elements: &sbdb_elements,
            }
            .position(tt, body, flags);
        }
        let (mut output, retflag) = if public_entry {
            local_jpl_position_with_public_entry(
                &leg,
                tt,
                original_ut.unwrap_or(tt),
                original_ut,
                body,
                flags,
                site,
                public_entry,
            )
        } else {
            local_jpl_position(
                &leg,
                tt,
                original_ut.unwrap_or(tt),
                original_ut,
                body,
                flags,
                site,
            )
        }?;
        local_jpl_sidereal_word(&mut output, tt, original_ut, flags, &|ut| {
            crate::planetary_moons::moon_spk_ayanamsha(ut, flags, &ayanamsha)
        })?;
        // Publish the acquired vector source, separately from analytical and
        // Keplerian exits. Explicit minor/satellite vectors use SPK; standard planets
        // retain their public source label. Exact-origin sentinels remain analytical.
        let source = if local_jpl_minor_body_candidate(body)
            || crate::planetary_moons::is_planetary_moon(body)
        {
            "SPK"
        } else if vector_self_observation(body, flags) || (body == SUN && flags & FLG_HELCTR != 0) {
            "Analytical"
        } else {
            "Skyfield"
        };
        native_source_result(Ok((output, retflag)), source)
    };
    let compute = |tt, body, flags| {
        let result =
            native_sidereal_frame_request(&compute_raw, &ayanamsha, selection, tt, body, flags);
        result.map_err(|error| directory.promote(original_ut.unwrap_or(tt), body, error))
    };
    let ut_to_tt = |ut| crate::time::default_timescale_tt_from_ut1(ut);
    let engine = CalcEngine {
        ut1_to_tt_jd: &ut_to_tt,
        compute_tt: &compute,
        nutation_obliquity_ut: &ecl_nut_result_ut,
        nutation_obliquity_tt: &|tt, flags| ecl_nut_result(tt, flags, &crate::erfa::nut06a),
        jpl_fallback_active: fallback_active,
    };
    let _trace_scope = NativeTraceInvocation::enter_factory();
    f(&engine, &|_, _, _| Ok(None))
}

/// Only a failed initial optional acquisition may select the legacy fallback.
/// A physical emission failure remains an error even when its category was
/// optional before the native central orbit existed.
fn optional_assist_result<T>(
    result: std::result::Result<Option<T>, crate::nbody::OwnedTrajectoryError>,
) -> Result<Option<T>> {
    use crate::nbody::OwnedTrajectoryError::{Emission, Initial};
    match result {
        Ok(value) => Ok(value),
        Err(Initial(error)) if crate::nbody::assist_fallback(&error) => Ok(None),
        Err(Initial(error) | Emission { error, .. }) => Err(error),
    }
}

/// Per-call model-acquisition variant retaining the corresponding static
/// provider/clock/flag pipeline. The borrowed callback is consumed synchronously
/// only at analytical geometric state reads; other bodies never call it.
#[allow(clippy::too_many_arguments)]
fn local_jpl_uranian_bary_state_with_models(
    models: Option<&AnalyticalModelCallback<'_>>,
    directory: &PublicModelDirectory<'_, '_>,
    tt: f64,
    body: i32,
) -> Result<PositionResult> {
    model_serving::uranian_bary_source_time_state_with_models(
        models,
        &|id, epoch| directory.state(id, epoch, 0.0),
        tt,
        body,
    )
}

/// One central support role's owned split-TDB record or typed-decline gap.
#[derive(Clone, Copy)]
struct PredictedSpkDomain {
    target: i32,
    anchor: crate::types::TwoPartEpoch,
    bounds: (f64, f64),
    relative: bool,
    declined: bool,
}

const PCTR_MAX_READ_ROLES: usize = 64;

/// Ordered central reads only; no physical state or coefficient is retained.
pub(in crate::positions) struct PredictedSpkReads {
    enabled: bool,
    domains: std::cell::RefCell<[Option<PredictedSpkDomain>; PCTR_MAX_READ_ROLES]>,
    role_limit: usize,
    count: std::cell::Cell<usize>,
    cursor: std::cell::Cell<usize>,
    recording: std::cell::Cell<bool>,
    complete: std::cell::Cell<bool>,
    active: std::cell::Cell<bool>,
    refused: std::cell::Cell<bool>,
    physical_error: std::cell::RefCell<Option<Error>>,
}

/// Keep the stencil's exact TT displacement, without reusing a rounded
/// central epoch or dividing a different represented spacing by the nominal h.
fn represented_predicted_sample_epoch(epoch: f64, offset: f64) -> Option<f64> {
    let at = epoch + offset;
    (offset == 0.0 || (at != epoch && at - epoch == offset)).then_some(at)
}

impl PredictedSpkReads {
    pub(in crate::positions) fn new(enabled: bool) -> Self {
        Self::with_role_limit(enabled, model_serving::PREDICTED_MAX_READ_ROLES)
    }

    /// PCTR has at most 12 state roles, each with three physical center probes
    /// and at most two stored-read routes. Fixed 64-slot stack metadata covers
    /// that 60-role recipe; ordinary predicted places retain their 32-role cap.
    pub(in crate::positions) fn pctr(enabled: bool) -> Self {
        Self::with_role_limit(enabled, PCTR_MAX_READ_ROLES)
    }

    fn with_role_limit(enabled: bool, role_limit: usize) -> Self {
        Self {
            enabled,
            domains: std::cell::RefCell::new([None; PCTR_MAX_READ_ROLES]),
            role_limit,
            count: std::cell::Cell::new(0),
            cursor: std::cell::Cell::new(0),
            recording: std::cell::Cell::new(true),
            complete: std::cell::Cell::new(true),
            active: std::cell::Cell::new(false),
            refused: std::cell::Cell::new(false),
            physical_error: std::cell::RefCell::new(None),
        }
    }

    /// Whether supplementary samples must retain their source epoch tails.
    pub(in crate::positions) fn is_refining(&self) -> bool {
        self.active.get()
    }

    /// Outward bounds for the difference of the two actual split epochs.
    /// These are TDB offsets, never a linearized TT-to-emission conversion.
    fn offsets(
        epoch: crate::types::TwoPartEpoch,
        anchor: crate::types::TwoPartEpoch,
    ) -> Option<(f64, f64)> {
        let whole = epoch.hi - anchor.hi;
        let fraction = epoch.lo - anchor.lo;
        if !whole.is_finite() || !fraction.is_finite() {
            return None;
        }
        let lower = (whole.next_down() + fraction.next_down()).next_down();
        let upper = (whole.next_up() + fraction.next_up()).next_up();
        (lower.is_finite() && upper.is_finite()).then_some((lower, upper))
    }

    #[cfg(test)]
    fn note(&self, target: i32, epoch: crate::types::TwoPartEpoch, offsets: (f64, f64)) {
        self.note_role(target, epoch, offsets, false, false);
    }

    fn note_role(
        &self,
        target: i32,
        epoch: crate::types::TwoPartEpoch,
        offsets: (f64, f64),
        relative: bool,
        declined: bool,
    ) {
        if !(offsets.0 < 0.0 && 0.0 < offsets.1) {
            self.complete.set(false);
            return;
        }
        let mut domains = self.domains.borrow_mut();
        let slot = self.count.get();
        if slot == self.role_limit {
            self.complete.set(false);
            return;
        }
        domains[slot] = Some(PredictedSpkDomain {
            target,
            anchor: epoch,
            bounds: offsets,
            relative,
            declined,
        });
        self.count.set(slot + 1);
    }

    fn role_admits(
        &self,
        role: usize,
        target: i32,
        epoch: crate::types::TwoPartEpoch,
        relative: bool,
    ) -> bool {
        self.domains
            .borrow()
            .get(role)
            .copied()
            .flatten()
            .filter(|domain| domain.target == target && domain.relative == relative)
            .and_then(|domain| {
                Self::offsets(epoch, domain.anchor)
                    .map(|delta| domain.bounds.0 < delta.0 && delta.1 < domain.bounds.1)
            })
            .unwrap_or(false)
    }

    /// Freeze immediately after the original central place succeeds. The
    /// legacy side places still execute, without changing these role records.
    pub(in crate::positions) fn finish_central(&self) {
        self.recording.set(false);
    }

    /// Data evaluation stays outside every bookkeeping borrow. Metadata is
    /// consulted only after an original successful read and cannot replace
    /// its result/error. Supplementary admission precedes the physical read.
    pub(in crate::positions) fn query<T>(
        &self,
        target: i32,
        epoch: crate::types::TwoPartEpoch,
        metadata: &dyn Fn(i32, f64, f64) -> Result<(f64, f64)>,
        read: impl FnOnce() -> Result<T>,
    ) -> Result<T> {
        self.query_owned(target, epoch, metadata, read, false)
    }

    /// Offset absence is a consumed source choice, bounded by its actual gap.
    pub(in crate::positions) fn relative_query<T>(
        &self,
        target: i32,
        epoch: crate::types::TwoPartEpoch,
        metadata: &dyn Fn(i32, f64, f64) -> Result<(f64, f64)>,
        read: impl FnOnce() -> Result<T>,
    ) -> Result<T> {
        self.query_owned(target, epoch, metadata, read, true)
    }

    /// Preserve an original offset-velocity probe without declaring its record
    /// a dependency of a position-only PPN place. Owned derivative samples do
    /// not acquire these unused channels. Original physical errors still make
    /// the central source unsuitable for supplementary differentiation.
    pub(in crate::positions) fn central_velocity_probe<T>(
        &self,
        epoch: crate::types::TwoPartEpoch,
        read: impl FnOnce() -> Result<T>,
    ) -> Result<T> {
        if self.active.get() {
            return self.refuse(epoch.value());
        }
        let value = read();
        if self.enabled
            && self.recording.get()
            && value
                .as_ref()
                .is_err_and(|error| !Self::optional_decline(error))
        {
            self.complete.set(false);
        }
        value
    }

    fn optional_decline(error: &Error) -> bool {
        matches!(
            error,
            Error::SpkOutOfRange { .. } | Error::SpkBodyNotFound { .. }
        )
    }

    /// A previous real read error wins over every later optional admission
    /// refusal. Once either has occurred, no further support I/O is permitted.
    pub(in crate::positions) fn before_read(&self, epoch: f64) -> Result<()> {
        if self.active.get() {
            if let Some(error) = self.physical_error.borrow().as_ref() {
                return Err(error.clone());
            }
            if self.refused.get() {
                return Err(Error::MathDomain { value: epoch });
            }
        }
        Ok(())
    }

    /// Internal admission marker, converted to a declined sample by `rates`.
    pub(in crate::positions) fn refuse<T>(&self, epoch: f64) -> Result<T> {
        self.refused.set(true);
        Err(Error::MathDomain { value: epoch })
    }

    /// Preserve real errors even if an optional PPN read handles them locally.
    pub(in crate::positions) fn physical_result<T>(&self, value: &Result<T>) {
        if self.active.get() && !self.refused.get() {
            if let Err(error) = value {
                if self.physical_error.borrow().is_none() {
                    *self.physical_error.borrow_mut() = Some(error.clone());
                }
            }
        }
    }

    fn query_owned<T>(
        &self,
        target: i32,
        epoch: crate::types::TwoPartEpoch,
        metadata: &dyn Fn(i32, f64, f64) -> Result<(f64, f64)>,
        read: impl FnOnce() -> Result<T>,
        optional_offset: bool,
    ) -> Result<T> {
        if !self.enabled {
            return read();
        }
        if self.active.get() {
            if let Some(error) = self.physical_error.borrow().as_ref() {
                return Err(error.clone());
            }
            if self.refused.get()
                || !self.role_admits(self.cursor.get(), target, epoch, optional_offset)
            {
                self.refused.set(true);
                return Err(Error::MathDomain {
                    value: epoch.value(),
                });
            }
            self.cursor.set(self.cursor.get() + 1);
        }
        let value = read();
        if self.active.get() {
            if let Err(error) = &value {
                if !(optional_offset && Self::optional_decline(error)) {
                    *self.physical_error.borrow_mut() = Some(error.clone());
                }
            }
            let declined = value
                .as_ref()
                .is_err_and(|error| optional_offset && Self::optional_decline(error));
            if self.domains.borrow()[self.cursor.get() - 1]
                .is_some_and(|domain| domain.declined != declined)
            {
                self.refused.set(true);
            }
        } else if self.recording.get()
            && value
                .as_ref()
                .is_err_and(|error| !(optional_offset && Self::optional_decline(error)))
        {
            // Even an old optional failure leaves its source uninferable.
            self.complete.set(false);
        } else if self.recording.get() && self.complete.get() {
            match metadata(target, epoch.hi, epoch.lo) {
                Ok(offsets) => {
                    self.note_role(target, epoch, offsets, optional_offset, value.is_err())
                }
                Err(_) => self.complete.set(false),
            }
        }
        value
    }

    fn rates(
        &self,
        place: &dyn Fn(f64) -> Result<[f64; 3]>,
        epoch: f64,
        outer_half_width: f64,
    ) -> Result<Option<[f64; 3]>> {
        self.rates_with_owned_leb(place, epoch, outer_half_width, false)
    }

    /// LEB-only support needs no fictitious SPK domain. Any new SPK target
    /// still fails matching central-role admission before reading; Native's
    /// nonempty gate is unchanged.
    pub(in crate::positions) fn rates_with_owned_leb(
        &self,
        place: &dyn Fn(f64) -> Result<[f64; 3]>,
        epoch: f64,
        outer_half_width: f64,
        owned_leb: bool,
    ) -> Result<Option<[f64; 3]>> {
        self.rates_with_direction_hint(
            place,
            epoch,
            outer_half_width,
            owned_leb,
            (f64::INFINITY, f64::INFINITY),
        )
    }

    /// Margins choose a candidate stencil direction only. They are not a TT
    /// interval inferred from TDB: every nonlinear read is still checked at
    /// its actual split epoch against the matching central role before I/O.
    pub(in crate::positions) fn rates_with_direction_hint(
        &self,
        place: &dyn Fn(f64) -> Result<[f64; 3]>,
        epoch: f64,
        outer_half_width: f64,
        owned_leb: bool,
        hint: (f64, f64),
    ) -> Result<Option<[f64; 3]>> {
        self.finish_central();
        if !self.enabled
            || !self.complete.get()
            || (!owned_leb && self.domains.borrow().iter().all(Option::is_none))
        {
            return Ok(None);
        }
        let (left, right) = self
            .domains
            .borrow()
            .iter()
            .flatten()
            .fold(hint, |(left, right), domain| {
                (left.min(-domain.bounds.0), right.min(domain.bounds.1))
            });
        let footprint = 2.0 * crate::orbital::NODAPS_SPEED_STEP_DAYS;
        let support = if left < footprint && left < right {
            (0.0, outer_half_width)
        } else if right < footprint && right < left {
            (-outer_half_width, 0.0)
        } else {
            (-outer_half_width, outer_half_width)
        };
        self.active.set(true);
        let result = model_serving::native_prediction_rates(
            &|offset| {
                self.refused.set(false);
                *self.physical_error.borrow_mut() = None;
                self.cursor.set(0);
                let Some(at) = represented_predicted_sample_epoch(epoch, offset) else {
                    self.refused.set(true);
                    return Ok(None);
                };
                let value = place(at);
                if let Some(error) = self.physical_error.borrow_mut().take() {
                    // A later role/count refusal cannot mask the first real
                    // source error, including one swallowed by optional PPN.
                    self.refused.set(false);
                    return Err(error);
                }
                if self.refused.get() || self.cursor.get() != self.count.get() {
                    self.refused.set(true);
                    return Ok(None);
                }
                value.map(Some)
            },
            support,
        );
        self.active.set(false);
        match result {
            Err(_) if self.refused.get() => Ok(None),
            _ => result.map(Some),
        }
    }
}

fn predicted_spk_rates_enabled(
    models: Option<&AnalyticalModelCallback<'_>>,
    body: i32,
    flags: i32,
) -> bool {
    models.is_none()
        && model_serving::is_builtin_predicted_rate_body(body)
        && flags & crate::constants::FLG_SPEED != 0
        && flags & (FLG_HELCTR | FLG_BARYCTR) == 0
}

#[allow(clippy::too_many_arguments)]
fn predicted_state_with_spk_rates(
    models: Option<&AnalyticalModelCallback<'_>>,
    state: &dyn Fn(i32, f64) -> Result<([f64; 3], [f64; 3])>,
    observation: &dyn Fn(f64) -> Result<model_serving::RetardedEarthObservation>,
    reads: &PredictedSpkReads,
    tt: f64,
    body: i32,
    flags: i32,
) -> Result<(PositionResult, model_serving::ModelRateOrigin)> {
    let rates =
        |place: &dyn Fn(f64) -> Result<[f64; 3]>, epoch, width| reads.rates(place, epoch, width);
    let central_complete = || reads.finish_central();
    model_serving::with_model_rate_origin(&rates, |rates| {
        model_serving::predicted_center_state_source_time_with_models_and_rates(
            models,
            state,
            observation,
            false,
            tt,
            body,
            flags,
            reads.enabled.then_some(rates),
            reads.enabled.then_some(&central_complete),
        )
    })
}

/// Per-call model-acquisition variant retaining the corresponding static
/// provider/clock/flag pipeline. The borrowed callback is consumed synchronously
/// only at analytical geometric state reads; other bodies never call it.
#[allow(clippy::too_many_arguments)]
fn local_jpl_predicted_center_state_with_models(
    models: Option<&AnalyticalModelCallback<'_>>,
    leg: &JplLeg<'_>,
    directory: &PublicModelDirectory<'_, '_>,
    tt: f64,
    body: i32,
    flags: i32,
) -> Result<(PositionResult, model_serving::ModelRateOrigin)> {
    let state = |id, epoch| {
        if id != SUN && id != EARTH && local_jpl_target(leg, id).is_err() {
            return Err(Error::UnknownBody { body_id: id });
        }
        directory.state(id, epoch, 0.0)
    };
    if !predicted_spk_rates_enabled(models, body, flags) {
        return model_serving::predicted_center_state_source_time_with_models(
            models,
            &state,
            &|epoch| directory.earth_observation(epoch),
            false,
            tt,
            body,
            flags,
        )
        .map(|state| (state, model_serving::ModelRateOrigin::LegacyStencil));
    }
    let reads = PredictedSpkReads::new(true);
    let metadata = |target, whole, fraction| {
        leg.registry()
            .geometric_ssb_record_offsets(target, whole, fraction)
    };
    let observed_state = |id, epoch| {
        let target = if id != SUN && id != EARTH {
            local_jpl_target(leg, id).map_err(|_| Error::UnknownBody { body_id: id })?
        } else {
            local_jpl_target(leg, id)?
        };
        // Repeat only the pure clock expression used by directory.state;
        // the physical read remains that original callback and split pair.
        let correction = crate::time::source_tdb_minus_tt_seconds(epoch, 0.0);
        let clock = crate::types::TwoPartEpoch {
            hi: epoch,
            lo: 0.0 + correction / 86400.0,
        };
        reads.query(target, clock, &metadata, || directory.state(id, epoch, 0.0))
    };
    let observed_earth = |epoch| {
        super::retarded_earth_observation_source_time_via(leg, epoch, &|target, whole, fraction| {
            reads.query(
                target,
                crate::types::TwoPartEpoch {
                    hi: whole,
                    lo: fraction,
                },
                &metadata,
                || directory.tdb(target, whole, fraction),
            )
        })
    };
    predicted_state_with_spk_rates(
        models,
        &observed_state,
        &observed_earth,
        &reads,
        tt,
        body,
        flags,
    )
}

/// Assemble the same callbacks with optional original public range metadata.
/// Legacy consumers supply None; only a genuine scalar entry installs a role.
pub(super) fn with_local_jpl_entry_outputs_and_fallback_and_range<R>(
    req: ServingRequest<'_>,
    original_ut: Option<f64>,
    fallback_active: bool,
    public_entry: bool,
    f: impl FnOnce(&CalcEngine<'_>, &CartesianCompute<'_>) -> Result<R>,
) -> Result<R> {
    with_local_jpl_entry_outputs_and_fallback_and_models(
        None,
        None,
        None,
        req,
        LocalJplEntry {
            original_ut,
            fallback_active,
            public_entry,
            owned_minor: None,
        },
        f,
    )
}

#[cfg(test)]
mod predicted_spk_rate_tests {
    use super::*;
    use crate::constants::{FLG_J2000, FLG_SPEED, FLG_TRUEPOS, NEPTUNE_LEVERRIER, PROSERPINA};
    use crate::types::TwoPartEpoch;
    use std::cell::{Cell, RefCell};

    #[test]
    fn pctr_dependency_capacity_does_not_expand_the_original_admission_limit() {
        for (reads, limit) in [
            (
                PredictedSpkReads::new(true),
                model_serving::PREDICTED_MAX_READ_ROLES,
            ),
            (PredictedSpkReads::pctr(true), PCTR_MAX_READ_ROLES),
        ] {
            for role in 0..limit {
                reads.note(
                    role as i32,
                    TwoPartEpoch {
                        hi: 2451545.0,
                        lo: 0.0,
                    },
                    (-1.0, 1.0),
                );
            }
            assert!(reads.complete.get());
            assert_eq!(reads.count.get(), limit);
            reads.note(
                999,
                TwoPartEpoch {
                    hi: 2451545.0,
                    lo: 0.0,
                },
                (-1.0, 1.0),
            );
            assert!(!reads.complete.get());
            assert_eq!(reads.count.get(), limit);
        }
    }

    #[test]
    fn emission_source_failure_never_becomes_an_optional_backend_miss() {
        let missing = || Error::EphemerisLoaderIo {
            path: "synthetic-force-file".into(),
            kind: std::io::ErrorKind::NotFound,
            errno: None,
            message: "missing".into(),
        };
        assert!(crate::nbody::assist_fallback(&missing()));
        assert_eq!(
            optional_assist_result::<()>(Err(crate::nbody::OwnedTrajectoryError::Initial(
                missing()
            )))
            .unwrap(),
            None
        );
        assert!(matches!(
            optional_assist_result::<()>(Err(
                crate::nbody::OwnedTrajectoryError::after_acquisition(missing())
            )),
            Err(Error::EphemerisLoaderIo { .. })
        ));
    }

    #[test]
    fn owned_observer_type3_tangent_uses_position_and_its_actual_clock() {
        let mut image = crate::spk::moon_test_kernel_bytes();
        image.resize(4096, 0);
        image[1048..1056].copy_from_slice(&(-172_800.0_f64).to_le_bytes());
        image[1056..1064].copy_from_slice(&172_800.0_f64.to_le_bytes());
        image[1064..1068].copy_from_slice(&399_i32.to_le_bytes());
        image[1068..1072].copy_from_slice(&0_i32.to_le_bytes());
        image[1072..1076].copy_from_slice(&1_i32.to_le_bytes());
        image[1076..1080].copy_from_slice(&3_i32.to_le_bytes());
        image[1080..1084].copy_from_slice(&385_i32.to_le_bytes());
        image[1084..1088].copy_from_slice(&416_i32.to_le_bytes());
        let au_km = 149_597_870.7;
        let slopes = [0.01, -0.02, 0.03];
        let mut words = Vec::new();
        for middle in [-86_400.0, 86_400.0] {
            words.extend([middle, 86_400.0]);
            for (axis, slope) in slopes.iter().enumerate() {
                words.extend([
                    (axis as f64 + 1.0 + slope * middle / 86400.0) * au_km,
                    slope * au_km,
                ]);
            }
            // Deliberately unrelated stored velocity channels, km/s.
            for _ in 0..3 {
                words.extend([99.0, 0.0]);
            }
        }
        words.extend([-172_800.0, 172_800.0, 14.0, 2.0]);
        for (index, value) in words.into_iter().enumerate() {
            let start = 3072 + index * 8;
            image[start..start + 8].copy_from_slice(&value.to_le_bytes());
        }
        let mut context =
            crate::context::EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
        let registry = context.spk_registry_mut();
        registry.register_image("observer-type3", image).unwrap();
        registry
            .register_body(EARTH, "observer-type3", 399)
            .unwrap();
        registry.register_body(399, "observer-type3", 399).unwrap();
        let provider = LocalSpkStateProvider::new(context.spk_registry(), context.mode()).unwrap();
        let central = provider
            .nodaps_state_two_part(
                EARTH,
                TwoPartEpoch {
                    hi: 2451545.0,
                    lo: 0.0,
                },
            )
            .unwrap();
        let derivative = crate::engine_api::NativeMinorAssist::earth_position_tt_rate(
            &context,
            2451545.0,
            central.pos(),
        )
        .unwrap();
        let clock_rate = crate::time::tdb_rate(2451545.0).unwrap();
        for axis in 0..3 {
            assert!(
                (derivative[axis] - slopes[axis] * clock_rate).abs() < slopes[axis].abs() * 1e-9,
                "{derivative:?}"
            );
            assert!((derivative[axis] - central.vel()[axis]).abs() > 0.001);
        }
    }

    #[test]
    fn declined_native_refinement_does_not_mark_its_eligible_state_owned() {
        let reads = PredictedSpkReads::new(true);
        reads.complete.set(false);
        let calls = RefCell::new(Vec::new());
        let state = |id, tt: f64| {
            calls.borrow_mut().push((id, tt.to_bits()));
            Ok(([0.0; 3], [0.0; 3]))
        };
        let observation = |tt: f64| {
            calls.borrow_mut().push((SUN, tt.to_bits()));
            Ok(model_serving::RetardedEarthObservation {
                relative: [1.0, 0.0, 0.0],
                sun_reception: [0.0; 3],
            })
        };
        let flags = FLG_SPEED | FLG_J2000 | FLG_TRUEPOS;
        let (actual, origin) = predicted_state_with_spk_rates(
            None,
            &state,
            &observation,
            &reads,
            2451545.0,
            NEPTUNE_LEVERRIER,
            flags,
        )
        .unwrap();
        let original_calls = calls.take();
        let original = model_serving::predicted_center_state_source_time_with_models(
            None,
            &state,
            &observation,
            false,
            2451545.0,
            NEPTUNE_LEVERRIER,
            flags,
        )
        .unwrap();
        assert_eq!(origin, model_serving::ModelRateOrigin::LegacyStencil);
        assert_eq!(actual.map(f64::to_bits), original.map(f64::to_bits));
        assert_eq!(original_calls, *calls.borrow());
    }

    fn epoch(offset: f64) -> TwoPartEpoch {
        TwoPartEpoch {
            hi: 2451545.0,
            lo: offset,
        }
    }

    #[test]
    fn central_read_roles_keep_distinct_records_of_the_same_target() {
        let reads = PredictedSpkReads::new(true);
        let read =
            |offset| reads.query(399, epoch(offset), &|_, _, _| Ok((-0.02, 0.02)), || Ok(()));
        // Reception and a retarded role consume disjoint real records. Each
        // role has its own smooth neighbourhood; their union is never a record.
        read(-0.125).unwrap();
        read(0.0).unwrap();
        let acquired = RefCell::new(Vec::new());
        let actual = reads
            .rates(
                &|tt| {
                    let dt = tt - 2451545.0;
                    for role in [-0.125, 0.0] {
                        reads.query(
                            399,
                            epoch(role + dt),
                            &|_, _, _| panic!("supplementary metadata is forbidden"),
                            || {
                                acquired.borrow_mut().push((role, dt));
                                Ok(())
                            },
                        )?;
                    }
                    Ok([120.0, 0.01 * (1.25 * dt).sin(), 40.0])
                },
                2451545.0,
                0.1,
            )
            .unwrap();
        let rates = actual.expect("distinct owned central roles remain differentiable");
        assert!((rates[1] - 0.0125).abs() <= 0.0125 * 1e-6);
        assert!(acquired.borrow().iter().all(|(_, dt)| dt.abs() < 0.02));
    }

    #[test]
    fn relative_offset_declines_keep_their_owned_gap_but_other_errors_do_not() {
        for error in [
            Error::SpkBodyNotFound { target: 599 },
            Error::SpkOutOfRange {
                target: 599,
                et: 0.0,
                start: 1.0,
                end: 2.0,
            },
        ] {
            let reads = PredictedSpkReads::new(true);
            let result: Result<()> =
                reads.relative_query(599, epoch(0.0), &|_, _, _| Ok((-0.02, 0.02)), || {
                    Err(error.clone())
                });
            assert_eq!(result.unwrap_err().to_string(), error.to_string());
            assert!(reads.complete.get());
            let acquired = RefCell::new(Vec::new());
            let rates = reads
                .rates_with_owned_leb(
                    &|tt| {
                        let dt = tt - 2451545.0;
                        let result: Result<()> = reads.relative_query(
                            599,
                            epoch(dt),
                            &|_, _, _| panic!("metadata is original-only"),
                            || {
                                acquired.borrow_mut().push(dt);
                                Err(error.clone())
                            },
                        );
                        // The normal center composer handles these exact two errors.
                        Ok(
                            if matches!(
                                result,
                                Err(Error::SpkOutOfRange { .. } | Error::SpkBodyNotFound { .. })
                            ) {
                                [120.0, 0.01 * (1.25 * dt).sin(), 40.0]
                            } else {
                                [120.0, 42.0, 40.0]
                            },
                        )
                    },
                    2451545.0,
                    0.1,
                    true,
                )
                .unwrap()
                .unwrap();
            assert!((rates[1] - 0.0125).abs() < 0.0125 * 1e-6);
            assert!(acquired.borrow().iter().all(|dt| dt.abs() < 0.02));
        }
        let reads = PredictedSpkReads::new(true);
        let result: Result<()> = reads.relative_query(
            599,
            epoch(0.0),
            &|_, _, _| panic!("failed data is not metadata"),
            || Err(Error::MathDomain { value: 123.0 }),
        );
        assert!(result.is_err());
        assert!(!reads.complete.get());
        assert_eq!(
            reads
                .rates_with_owned_leb(&|_| panic!("no sample"), 2451545.0, 0.1, true)
                .unwrap(),
            None
        );
    }

    #[test]
    fn unused_central_velocity_probes_preserve_errors_without_source_roles() {
        for error in [
            Error::SpkBodyNotFound { target: 599 },
            Error::MathDomain { value: 123.0 },
        ] {
            let reads = PredictedSpkReads::new(true);
            let value: Result<()> = reads.central_velocity_probe(epoch(0.0), || Err(error.clone()));
            assert_eq!(value.unwrap_err(), error);
            assert_eq!(reads.count.get(), 0);
            assert_eq!(
                reads.complete.get(),
                matches!(error, Error::SpkBodyNotFound { .. })
            );
        }
    }

    #[test]
    fn earlier_physical_error_wins_over_later_leb_admission_and_stops_reads() {
        for physical_first in [true, false] {
            let reads = PredictedSpkReads::new(true);
            let physical_reads = Cell::new(0);
            let result = reads.rates_with_owned_leb(
                &|tt| {
                    if physical_first {
                        reads.before_read(tt)?;
                        physical_reads.set(physical_reads.get() + 1);
                        reads.physical_result(&Err::<(), _>(Error::MathDomain { value: 777.0 }));
                        // The optional deflector swallowed its read error. A later
                        // scalar call must check this guard BEFORE refusing LEB.
                        let later = reads.before_read(tt + 100.0);
                        if later.is_ok() {
                            let _ = reads.refuse::<()>(tt + 100.0);
                        }
                        assert!(matches!(later, Err(Error::MathDomain { value: 777.0 })));
                    } else {
                        let _ = reads.refuse::<()>(tt);
                        if reads.before_read(tt + 100.0).is_ok() {
                            physical_reads.set(physical_reads.get() + 1);
                            reads
                                .physical_result(&Err::<(), _>(Error::MathDomain { value: 777.0 }));
                        }
                    }
                    Ok([120.0, 0.0, 40.0])
                },
                2451545.0,
                0.1,
                true,
            );
            if physical_first {
                assert!(matches!(result, Err(Error::MathDomain { value: 777.0 })));
                assert_eq!(physical_reads.get(), 1);
            } else {
                assert_eq!(result.unwrap(), None);
                assert_eq!(physical_reads.get(), 0);
            }
        }
    }

    #[test]
    fn split_domains_intersect_per_actual_target_and_refuse_extra_targets() {
        let reads = PredictedSpkReads::new(true);
        reads.note(599, epoch(0.0), (-0.5, 0.5));
        reads.note(599, epoch(0.125), (-0.25, 0.5));
        assert!(reads.role_admits(0, 599, epoch(0.0), false));
        assert!(reads.role_admits(1, 599, epoch(0.25), false));
        // Same target, distinct roles: do not intersect their records.
        assert!(reads.role_admits(0, 599, epoch(-0.125), false));
        assert!(!reads.role_admits(1, 599, epoch(-0.125), false));
        assert!(!reads.role_admits(0, 599, epoch(0.5), false));
        assert!(!reads.role_admits(0, 5, epoch(0.0), false));
        for target in [10, 399, 699] {
            reads.note(target, epoch(0.0), (-1.0, 1.0));
        }
        assert!(reads.complete.get());
        assert!(!reads.role_admits(5, 6, epoch(0.0), false));
        assert_eq!(reads.domains.borrow().iter().flatten().count(), 5);
        for _ in reads.count.get()..model_serving::PREDICTED_MAX_READ_ROLES {
            reads.note(10, epoch(0.0), (-1.0, 1.0));
        }
        assert!(reads.complete.get());
        reads.note(10, epoch(0.0), (-1.0, 1.0));
        assert!(!reads.complete.get());
        assert_eq!(
            reads.domains.borrow().iter().flatten().count(),
            model_serving::PREDICTED_MAX_READ_ROLES
        );
        let low = 2.0_f64.powi(-40);
        assert_eq!(epoch(low).value(), epoch(0.0).value());
        let narrow = PredictedSpkReads::new(true);
        narrow.note(399, epoch(0.0), (-low / 2.0, low / 2.0));
        assert!(!narrow.role_admits(0, 399, epoch(low), false));
    }

    #[test]
    fn joins_and_metadata_errors_retain_the_original_result_without_side_reads() {
        for offsets in [Ok((0.0, 1.0)), Err(Error::MathDomain { value: 123.0 })] {
            let reads = PredictedSpkReads::new(true);
            let expected = ([1.0, -0.0, 3.0], [4.0, 5.0, 6.0]);
            assert_eq!(
                reads
                    .query(399, epoch(0.0), &|_, _, _| offsets.clone(), || Ok(expected))
                    .unwrap(),
                expected
            );
            assert!(!reads.complete.get());
            assert_eq!(
                reads
                    .rates(
                        &|_| panic!("uninferable domain must not sample"),
                        2451545.0,
                        0.1
                    )
                    .unwrap(),
                None
            );
        }
    }

    #[test]
    fn swallowed_optional_admission_refusal_never_evaluates_an_unowned_state() {
        let reads = PredictedSpkReads::new(true);
        reads.note(10, epoch(0.0), (-0.02, 0.02));
        let acquired = RefCell::new(Vec::new());
        let rates = reads
            .rates(
                &|tt| {
                    let offset = tt - 2451545.0;
                    let result = reads.query(
                        10,
                        epoch(offset),
                        &|_, _, _| panic!("no supplementary metadata"),
                        || {
                            acquired.borrow_mut().push(offset);
                            Ok(([0.0; 3], [0.0; 3]))
                        },
                    );
                    // This models the pre-existing optional-deflector error policy.
                    // Its omission value must not enter an admitted derivative.
                    Ok(if result.is_err() {
                        [120.0, 42.0, 40.0]
                    } else {
                        [120.0, 0.01 * (1.25 * offset).sin(), 40.0]
                    })
                },
                2451545.0,
                0.1,
            )
            .unwrap()
            .unwrap();
        assert!((rates[1] - 0.0125).abs() < 0.0125 * 1e-6);
        assert!(acquired.borrow().iter().all(|offset| offset.abs() < 0.02));
    }

    #[test]
    fn swallowed_supplementary_physical_error_propagates_without_more_data_reads() {
        let reads = PredictedSpkReads::new(true);
        reads.note(10, epoch(0.0), (-1.0, 1.0));
        let acquired = Cell::new(0);
        let result = reads.rates(
            &|tt| {
                let source = || -> Result<([f64; 3], [f64; 3])> {
                    acquired.set(acquired.get() + 1);
                    Err(Error::MathDomain { value: 777.0 })
                };
                let _ = reads.query(10, epoch(tt - 2451545.0), &|_, _, _| unreachable!(), source);
                let _ = reads.query(10, epoch(tt - 2451545.0), &|_, _, _| unreachable!(), source);
                Ok([120.0, 1.0, 40.0])
            },
            2451545.0,
            0.1,
        );
        assert_eq!(result, Err(Error::MathDomain { value: 777.0 }));
        assert_eq!(acquired.get(), 1);
    }

    fn closed_state(
        enabled: bool,
        models: Option<&AnalyticalModelCallback<'_>>,
        body: i32,
    ) -> (PositionResult, Vec<(i32, u64, u64)>, usize) {
        closed_state_on_records(enabled, models, body, false)
    }

    fn closed_state_on_records(
        enabled: bool,
        models: Option<&AnalyticalModelCallback<'_>>,
        body: i32,
        separate: bool,
    ) -> (PositionResult, Vec<(i32, u64, u64)>, usize) {
        let reads = PredictedSpkReads::new(enabled);
        let acquired = RefCell::new(Vec::new());
        let metadata_count = Cell::new(0);
        let metadata = |_: i32, whole: f64, fraction: f64| {
            metadata_count.set(metadata_count.get() + 1);
            let (before, after) = if separate {
                let start =
                    (((whole - 2451545.0) + fraction + 0.015625) * 16.0).floor() / 16.0 - 0.015625;
                (2451545.0 + start, 2451545.0 + start + 0.0625)
            } else {
                (2451543.0, 2451547.0)
            };
            Ok(((before - whole) - fraction, (after - whole) - fraction))
        };
        // A closed synthetic clock F(TT)=TT+1/8 and explicit polynomial
        // domains certify these invented sources, not opaque public callbacks.
        let at = |target, tt: f64, value| {
            reads.query(
                target,
                TwoPartEpoch { hi: tt, lo: 0.125 },
                &metadata,
                || {
                    acquired
                        .borrow_mut()
                        .push((target, tt.to_bits(), 0.125_f64.to_bits()));
                    Ok((value, [0.0; 3]))
                },
            )
        };
        let state = |body, tt| {
            assert_eq!(body, EARTH);
            at(399, tt, [0.0; 3])
        };
        let observation = |tt: f64| {
            let sun = at(10, tt, [0.0; 3])?.0;
            let model = if body == PROSERPINA {
                crate::hypothetical::calc_proserpina(tt)?
            } else {
                crate::hypothetical::calc_fictitious_position(body, tt)?
            };
            let to_j2000 = |lon, lat, radius| -> Result<[f64; 3]> {
                let (lon, lat) = if body == PROSERPINA {
                    crate::hypothetical::precess_ecliptic(lon, lat, tt, 2451545.0)?
                } else {
                    (lon, lat)
                };
                Ok(model_serving::model_xyz(lon, lat, radius))
            };
            let target = to_j2000(model.lon_deg, model.lat_deg, model.dist_au)?;
            let offset = tt - 2451545.0;
            let desired = to_j2000(
                120.0 + 0.02 * (0.5 * offset).sin(),
                0.01 * (1.25 * offset).sin(),
                40.0 + 0.01 * (0.75 * offset).sin(),
            )?;
            let relative: [f64; 3] = std::array::from_fn(|i| target[i] - desired[i]);
            let (sine, cosine) = (84381.448_f64 / 3600.0).to_radians().sin_cos();
            let earth = at(
                399,
                tt,
                [
                    relative[0],
                    cosine * relative[1] - sine * relative[2],
                    sine * relative[1] + cosine * relative[2],
                ],
            )?
            .0;
            Ok(model_serving::RetardedEarthObservation {
                relative: earth,
                sun_reception: sun,
            })
        };
        let (result, rate_origin) = predicted_state_with_spk_rates(
            models,
            &state,
            &observation,
            &reads,
            2451545.0,
            body,
            FLG_J2000 | FLG_SPEED | FLG_TRUEPOS,
        )
        .unwrap();
        assert_eq!(
            rate_origin,
            if enabled {
                model_serving::ModelRateOrigin::OwnedDerivative
            } else {
                model_serving::ModelRateOrigin::LegacyStencil
            }
        );
        (result, acquired.into_inner(), metadata_count.get())
    }

    #[test]
    fn owned_spk_model_differentiates_all_three_closed_components_after_original_reads() {
        let original = closed_state(false, None, NEPTUNE_LEVERRIER);
        let refined = closed_state(true, None, NEPTUNE_LEVERRIER);
        assert_eq!(
            refined.0[..3]
                .iter()
                .map(|x| x.to_bits())
                .collect::<Vec<_>>(),
            original.0[..3]
                .iter()
                .map(|x| x.to_bits())
                .collect::<Vec<_>>()
        );
        assert!(refined.1.starts_with(&original.1));
        assert_eq!(original.1.len(), 9);
        assert!(refined.1.len() > original.1.len());
        assert_eq!(refined.2, 3); // Only the three central roles own metadata.
        for (actual, expected) in refined.0[3..].iter().zip([0.01_f64, 0.0125, 0.0075]) {
            assert!(
                (actual - expected).abs() <= expected.abs() * 1e-6,
                "actual={actual}, closed={expected}"
            );
        }
    }

    #[test]
    fn central_roles_refine_after_legacy_side_places_cross_record_joins() {
        let original = closed_state_on_records(false, None, NEPTUNE_LEVERRIER, true);
        let refined = closed_state_on_records(true, None, NEPTUNE_LEVERRIER, true);
        assert_eq!(refined.0[..3], original.0[..3]);
        assert!(refined.1.starts_with(&original.1));
        assert_eq!(original.1.len(), 9);
        assert_eq!(refined.2, 3);
        for (actual, expected) in refined.0[3..].iter().zip([0.01_f64, 0.0125, 0.0075]) {
            assert!((actual - expected).abs() <= expected.abs() * 1e-6);
        }
    }

    #[test]
    fn narrow_owned_side_chooses_one_sided_stencil_but_checks_nonlinear_clock() {
        let reads = PredictedSpkReads::new(true);
        reads.note(399, epoch(0.0), (-1e-9, 0.04));
        let acquired = RefCell::new(Vec::new());
        let actual = reads
            .rates(
                &|tt| {
                    let dt = tt - 2451545.0;
                    let source_offset = dt + 0.25 * dt * dt;
                    reads.query(399, epoch(source_offset), &|_, _, _| unreachable!(), || {
                        acquired.borrow_mut().push(source_offset);
                        Ok(())
                    })?;
                    Ok([120.0, 0.01 * (1.25 * dt).sin(), 40.0])
                },
                2451545.0,
                0.1,
            )
            .unwrap()
            .unwrap();
        assert!((actual[1] - 0.0125).abs() <= 0.0125 * 1e-6);
        assert!(acquired.borrow().contains(&0.0));
        assert!(
            acquired
                .borrow()
                .iter()
                .all(|&offset| (0.0..0.04).contains(&offset))
        );
    }

    #[test]
    fn sub_ulp_owned_window_declines_without_repeating_scalar_tt_reads() {
        let reads = PredictedSpkReads::new(true);
        let physical_reads = Cell::new(0);
        let half_width = (2451545.0_f64.next_up() - 2451545.0) / 8.0;
        reads
            .query(
                399,
                epoch(0.0),
                &|_, _, _| Ok((-half_width, half_width)),
                || {
                    physical_reads.set(physical_reads.get() + 1);
                    Ok(())
                },
            )
            .unwrap();
        let result = reads
            .rates(
                &|at| {
                    let dt = at - 2451545.0;
                    reads.query(399, epoch(dt), &|_, _, _| unreachable!(), || {
                        physical_reads.set(physical_reads.get() + 1);
                        Ok(())
                    })?;
                    // The underlying closed rate is 0.25 deg/day, not zero. Its
                    // neighbouring scalar TT values lie outside this source window.
                    Ok([120.0, 0.25 * dt, 40.0])
                },
                2451545.0,
                0.1,
            )
            .unwrap();
        assert_eq!(result, None);
        assert_eq!(
            physical_reads.get(),
            1,
            "only the original central read is allowed"
        );
    }

    #[test]
    fn represented_sample_clock_rejects_collapsed_and_changed_nonzero_spacing() {
        let epoch = 2451545.0_f64;
        let ulp = epoch.next_up() - epoch;
        assert_eq!(represented_predicted_sample_epoch(epoch, 0.0), Some(epoch));
        assert_eq!(
            represented_predicted_sample_epoch(epoch, ulp),
            Some(epoch.next_up())
        );
        assert_eq!(
            represented_predicted_sample_epoch(epoch, -ulp),
            Some(epoch.next_down())
        );
        assert_eq!(represented_predicted_sample_epoch(epoch, ulp / 4.0), None);
        let changed = ulp * 0.75;
        assert_ne!(epoch + changed, epoch);
        assert_ne!((epoch + changed) - epoch, changed);
        assert_eq!(represented_predicted_sample_epoch(epoch, changed), None);
    }

    #[test]
    fn role_order_count_kind_and_decline_status_cannot_change_during_refinement() {
        for mutation in ["order", "missing", "extra", "kind", "decline"] {
            let reads = PredictedSpkReads::new(true);
            reads
                .query(399, epoch(0.0), &|_, _, _| Ok((-1.0, 1.0)), || Ok(()))
                .unwrap();
            reads
                .relative_query(599, epoch(0.0), &|_, _, _| Ok((-1.0, 1.0)), || Ok(()))
                .unwrap();
            let calls = Cell::new(0);
            let source = || {
                calls.set(calls.get() + 1);
                Ok(())
            };
            let actual = reads
                .rates(
                    &|tt| {
                        let before = calls.get();
                        let at = epoch(tt - 2451545.0);
                        if mutation == "order" {
                            let _ =
                                reads.relative_query(599, at, &|_, _, _| unreachable!(), source);
                        } else if mutation == "kind" {
                            let _ =
                                reads.relative_query(399, at, &|_, _, _| unreachable!(), source);
                        } else {
                            reads.query(399, at, &|_, _, _| unreachable!(), source)?;
                            if mutation != "missing" {
                                let _ = reads.relative_query(
                                    599,
                                    at,
                                    &|_, _, _| unreachable!(),
                                    || {
                                        calls.set(calls.get() + 1);
                                        if mutation == "decline" {
                                            Err(Error::SpkBodyNotFound { target: 599 })
                                        } else {
                                            Ok(())
                                        }
                                    },
                                );
                            }
                            if mutation == "extra" {
                                let _ = reads.query(399, at, &|_, _, _| unreachable!(), source);
                            }
                        }
                        assert!(calls.get() - before <= 2, "mutation={mutation}");
                        Ok([120.0, 1.0, 40.0])
                    },
                    2451545.0,
                    0.1,
                )
                .unwrap();
            assert_eq!(actual, None, "mutation={mutation}");
            if matches!(mutation, "order" | "kind") {
                assert_eq!(calls.get(), 0);
            }
        }
    }

    #[test]
    fn owned_proserpina_rates_resolve_the_closed_moving_frame_derivative() {
        // The observer leaves a prescribed spherical motion in mean date axes.
        // Its closed derivative tests the complete moving-frame place, rather
        // than treating this date-defined model as a fixed J2000 orbit.
        assert!(predicted_spk_rates_enabled(None, PROSERPINA, FLG_SPEED));
        let original = closed_state(false, None, PROSERPINA);
        let refined = closed_state(true, None, PROSERPINA);
        assert_eq!(refined.0[..3], original.0[..3]);
        assert!(refined.1.starts_with(&original.1));
        assert!(refined.1.len() > original.1.len());
        for (axis, expected) in [0.01_f64, 0.0125, 0.0075].into_iter().enumerate() {
            assert!((original.0[axis + 3] - expected).abs() > expected.abs() * 1e-6);
            assert!(
                (refined.0[axis + 3] - expected).abs() <= expected.abs() * 1e-6,
                "component {axis}: actual={}, exact={expected}",
                refined.0[axis + 3],
            );
        }
    }

    #[test]
    fn opaque_callbacks_do_not_acquire_owned_domains() {
        let model = |_, _, _| Ok(None);
        assert!(!predicted_spk_rates_enabled(
            Some(&model),
            NEPTUNE_LEVERRIER,
            FLG_SPEED
        ));
        assert!(!predicted_spk_rates_enabled(
            Some(&model),
            PROSERPINA,
            FLG_SPEED
        ));
        for flags in [0, FLG_SPEED | FLG_HELCTR, FLG_SPEED | FLG_BARYCTR] {
            assert!(!predicted_spk_rates_enabled(None, NEPTUNE_LEVERRIER, flags));
        }
        let enabled = predicted_spk_rates_enabled(Some(&model), NEPTUNE_LEVERRIER, FLG_SPEED);
        let callback = closed_state(enabled, Some(&model), NEPTUNE_LEVERRIER);
        let legacy = closed_state(false, Some(&model), NEPTUNE_LEVERRIER);
        assert_eq!(callback, legacy);
        assert_eq!(callback.2, 0);
        for body in [crate::constants::VULCAN, PROSERPINA] {
            let model_calls = RefCell::new(Vec::new());
            let support_calls = RefCell::new(Vec::new());
            let model = |family, tt: f64, source_time| {
                model_calls
                    .borrow_mut()
                    .push((family, tt.to_bits(), source_time));
                Ok(Some(crate::hypothetical::HypoState {
                    lon_deg: 45.0,
                    lat_deg: 3.0,
                    dist_au: 2.0,
                    dlon: 0.0,
                    dlat: 0.0,
                    ddist: 0.0,
                }))
            };
            let models = Some(&model as &AnalyticalModelCallback<'_>);
            let flags = FLG_J2000 | FLG_SPEED | FLG_TRUEPOS;
            let reads = PredictedSpkReads::new(predicted_spk_rates_enabled(models, body, flags));
            assert!(!reads.enabled);
            let state = |id, tt: f64| {
                support_calls.borrow_mut().push((id, tt.to_bits()));
                Ok(([0.0; 3], [0.0; 3]))
            };
            let observation = |tt: f64| {
                support_calls.borrow_mut().push((SUN, tt.to_bits()));
                Ok(model_serving::RetardedEarthObservation {
                    relative: [1.0, 0.0, 0.0],
                    sun_reception: [0.0; 3],
                })
            };
            let actual = predicted_state_with_spk_rates(
                models,
                &state,
                &observation,
                &reads,
                2451545.0,
                body,
                flags,
            )
            .unwrap();
            let actual_models = model_calls.take();
            let actual_support = support_calls.take();
            let original = model_serving::predicted_center_state_source_time_with_models(
                models,
                &state,
                &observation,
                false,
                2451545.0,
                body,
                flags,
            )
            .unwrap();
            assert_eq!(actual.0.map(f64::to_bits), original.map(f64::to_bits));
            assert_eq!(actual.1, model_serving::ModelRateOrigin::LegacyStencil);
            assert_eq!(actual_models, *model_calls.borrow());
            assert_eq!(actual_support, *support_calls.borrow());
            assert_eq!(actual_models.is_empty(), body == PROSERPINA);
            assert!(reads.domains.borrow().iter().all(Option::is_none));
        }
    }
}
