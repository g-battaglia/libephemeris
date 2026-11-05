//! Planet-centered positions and rates over explicitly acquired endpoint states.
//!
//! Scalar TT and retained TT/TDB epochs feed the same three-step target
//! retardation, enabled deflection and aberration, and output-frame reduction.
//! Spherical positions use degrees and AU; rates are per day. Output flags
//! select Cartesian components or radians at the existing final reduction.
//! Provider, coverage, numerical and sidereal-configuration errors propagate;
//! this module does not select replacement data or alter acquisition policy.

use super::{
    echo_request_bits, exclusive_ephemeris_bit, finalize_output_flags, pctr_ayanamsa_deg,
    strip_output_flags,
};
use crate::constants::{
    FLG_BARYCTR, FLG_EQUATORIAL, FLG_HELCTR, FLG_ICRS, FLG_J2000, FLG_NOABERR, FLG_NOGDEFL,
    FLG_NONUT, FLG_SIDEREAL, FLG_SPEED, FLG_SPEED3, FLG_TRUEPOS, FLG_XYZ, MOON,
};
use crate::context::{EphemerisContext, ProviderMode};
use crate::errors::{Error, Result};
use crate::positions::routing::PredictedSpkReads;
use crate::positions::{
    BODY_SPEED_HALF_STEP_DAYS, C_AU_PER_DAY_PCTR, GeometricState, MOON_SPEED_HALF_STEP_DAYS,
    PctrEndpoint, PctrEngine, PositionResult, StateFrame, StateProvider, apply_aberration,
    apply_gravitational_deflection_with_lookup, cartesian_to_spherical_deg, check_finite_vec3,
    context_serving, frame_matrices_for_date, frame_transform, mat3_vec3, nutation_rate_deg_day,
    reference_frames, rotate_equatorial_to_ecliptic, vec3_dist, vondrak_mean_obliquity_rad,
    vondrak_precession_matrix,
};
use crate::types::TwoPartEpoch;
use std::cell::RefCell;

#[cfg(test)]
use context_serving::hypothetical_support_observations as observations;

/// Select the initial planet-centric rate half-step: six seconds for the Moon, one otherwise.
pub fn pctr_speed_half_step_days(body_id: i32) -> f64 {
    if body_id == MOON {
        MOON_SPEED_HALF_STEP_DAYS
    } else {
        BODY_SPEED_HALF_STEP_DAYS
    }
}

/// Whether a body is a supported planet-centric endpoint: Sun, Moon,
/// Mercury through Pluto or Earth. Other model families are not endpoints.
pub fn is_pctr_supported_body(body_id: i32) -> bool {
    super::body_uses_jpl_ephemeris(body_id)
}

/// Planet-centric apparent ICRS ray in AU and one-way light time in days.
///
/// The observer stays at reception while the target undergoes three light-time
/// iterations. TRUEPOS returns instantaneous geometry instead. Enabled
/// planetary deflection then aberration supply the apparent ray. Bound-state
/// acquisition and correction errors propagate.
pub fn pctr_apparent_vector(
    engine: &PctrEngine<'_>,
    jd_tt: f64,
    calc_flags: i32,
) -> Result<([f64; 3], f64)> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let (tgt_p, tgt_v) = (engine.state_at)(PctrEndpoint::Target, jd_tt)?;
    let (obs_p, obs_v) = (engine.state_at)(PctrEndpoint::Observer, jd_tt)?;
    check_finite_vec3(tgt_p)?;
    check_finite_vec3(tgt_v)?;
    check_finite_vec3(obs_p)?;
    check_finite_vec3(obs_v)?;
    let mut p = [
        tgt_p[0] - obs_p[0],
        tgt_p[1] - obs_p[1],
        tgt_p[2] - obs_p[2],
    ];
    if (calc_flags & FLG_TRUEPOS) != 0 {
        return Ok((p, 0.0));
    }
    let mut light_time = 0.0;
    for _ in 0..3 {
        let dist_au = vec3_dist(p);
        light_time = dist_au / C_AU_PER_DAY_PCTR;
        let (ret_p, ret_v) = (engine.state_at)(PctrEndpoint::Target, jd_tt - light_time)?;
        check_finite_vec3(ret_p)?;
        check_finite_vec3(ret_v)?;
        p = [
            ret_p[0] - obs_p[0],
            ret_p[1] - obs_p[1],
            ret_p[2] - obs_p[2],
        ];
    }
    if (calc_flags & FLG_NOGDEFL) == 0 {
        p = (engine.apply_deflection)(p, obs_p, jd_tt, light_time)?;
    }
    if (calc_flags & FLG_NOABERR) == 0 {
        p = (engine.apply_aberration)(p, obs_v, light_time)?;
    }
    Ok((p, light_time))
}

/// Project a planet-centric apparent ray through its output frame.
/// Return longitude/latitude in degrees and radius in AU, with the selected
/// ecliptic sidereal offset. Provider, frame and anchor errors propagate.
pub fn pctr_position_no_speed(
    engine: &PctrEngine<'_>,
    jd_tt: f64,
    calc_flags: i32,
) -> Result<(f64, f64, f64)> {
    let (p, _) = pctr_apparent_vector(engine, jd_tt, calc_flags)?;
    let (mut p1, p2, p3) = (engine.frame_at)(p, jd_tt, calc_flags)?;
    if !p1.is_finite() || !p2.is_finite() || !p3.is_finite() {
        return Err(Error::NonFiniteArgument { value: p1 });
    }
    if (calc_flags & FLG_SIDEREAL) != 0 && (calc_flags & FLG_EQUATORIAL) == 0 {
        let ayanamsha = pctr_ayanamsa_deg(engine, jd_tt)?;
        p1 = (p1 - ayanamsha).rem_euclid(360.0);
    }
    Ok((p1, p2, p3))
}

/// Differentiate the complete planet-centric no-speed position at symmetric
/// TT offsets. Unwrap longitude over its short arc; ecliptic sidereal rates
/// subtract the separately sampled ayanamsha drift. Equatorial sidereal rates
/// differentiate their mean-equator output directly. Errors propagate.
pub fn pctr_central_speed(
    engine: &PctrEngine<'_>,
    jd_tt: f64,
    body_id: i32,
    calc_flags: i32,
    equatorial: bool,
) -> Result<(f64, f64, f64)> {
    let dt = pctr_speed_half_step_days(body_id);
    let mut sample_flags = calc_flags & !(FLG_SPEED | FLG_SPEED3);
    if !equatorial {
        sample_flags &= !FLG_SIDEREAL;
    }
    let (p1_prev, p2_prev, p3_prev) = pctr_position_no_speed(engine, jd_tt - dt, sample_flags)?;
    let (p1_next, p2_next, p3_next) = pctr_position_no_speed(engine, jd_tt + dt, sample_flags)?;
    let step = 2.0 * dt;
    let mut dp1 = (p1_next - p1_prev) / step;
    let wrap_threshold = 180.0 / step;
    if dp1 > wrap_threshold {
        dp1 -= 360.0 / step;
    } else if dp1 < -wrap_threshold {
        dp1 += 360.0 / step;
    }
    let dp2 = (p2_next - p2_prev) / step;
    let dp3 = (p3_next - p3_prev) / step;
    if (calc_flags & FLG_SIDEREAL) != 0 && !equatorial {
        let ayanamsha_prev = pctr_ayanamsa_deg(engine, jd_tt - dt)?;
        let ayanamsha_next = pctr_ayanamsa_deg(engine, jd_tt + dt)?;
        let delta = (ayanamsha_next - ayanamsha_prev + 180.0).rem_euclid(360.0) - 180.0;
        dp1 -= delta / step;
    }
    Ok((dp1, dp2, dp3))
}

/// Calculate a standard planet as observed from another standard planet.
///
/// `tjdet` is TT JD and distance is target-to-observer AU. Unsupported body
/// identifiers return [`Error::UnknownBody`]. HELCTR/BARYCTR request an
/// astrometric place without changing the explicit observer; TOPOCTR is inert.
/// The result uses the selected public representation and returned flags.
/// Provider, numerical, frame and sidereal-model errors propagate.
pub fn calc_pctr(
    engine: &PctrEngine<'_>,
    tjdet: f64,
    target: i32,
    center: i32,
    flags: i32,
) -> Result<(PositionResult, i32)> {
    if !tjdet.is_finite() {
        return Err(Error::NonFiniteArgument { value: tjdet });
    }
    let mut mapped = flags;
    if (mapped & (FLG_HELCTR | FLG_BARYCTR)) != 0 {
        mapped = (mapped | FLG_NOABERR | FLG_NOGDEFL) & !(FLG_HELCTR | FLG_BARYCTR);
    }
    if !is_pctr_supported_body(target) {
        return Err(Error::UnknownBody { body_id: target });
    }
    if !is_pctr_supported_body(center) {
        return Err(Error::UnknownBody { body_id: center });
    }
    let norm_flags = exclusive_ephemeris_bit(mapped, engine.jpl_fallback_active);
    let calc_iflag = strip_output_flags(norm_flags);
    let equatorial = (calc_iflag & FLG_EQUATORIAL) != 0;
    let (p1, p2, p3) = pctr_position_no_speed(engine, tjdet, calc_iflag)?;
    // The pctr convention gates velocity on FLG_SPEED only: SPEED3 alone
    // is echoed but leaves the velocity slots zero.
    let (dp1, dp2, dp3) = if (calc_iflag & FLG_SPEED) != 0 {
        pctr_central_speed(engine, tjdet, target, calc_iflag, equatorial)?
    } else {
        (0.0, 0.0, 0.0)
    };
    Ok(finalize_output_flags(
        [p1, p2, p3, dp1, dp2, dp3],
        calc_iflag,
        norm_flags,
    ))
}

/// Context-configured planet-centric position at a finite scalar TT Julian Day.
///
/// Standard target and observer IDs are Sun through Pluto plus Earth. HELCTR
/// and BARYCTR select astrometric corrections while retaining the requested
/// observer; TOPOCTR has no surface-observer effect. Position units are
/// degrees/degrees/AU and rates per day, with XYZ/RADIANS applied once.
///
/// Explicit NativeJpl uses the context's owned local SPK registry and split
/// TT/TDB sampling. SealedLeb/Auto use their bound LEB stack, whose vector
/// reader consumes scalar TT independently at every read. No network or
/// provider substitution occurs. Three target-only retardation iterations
/// retain the observation TDB split; SPEED retains one-/six-second TT sample
/// fractions. SPEED3 alone returns zero rates. The source's native vector
/// nutation and independent default-timescale UT1 supply the frame and sidereal
/// clocks. Caller Delta-T selects sidereal defining epochs, not observation TT.
/// Uncomposed ICRS LEB endpoints supply stored TT tangents for geometric rates,
/// including stored system barycenters without a planet-center correction.
/// Native Type 2/3 position-polynomial chains supply separate TDB tangents with
/// the split source clock derivative, strictly inside record/selection domains.
/// Stored Type 3 physical velocity is not used as its position derivative.
/// Other geometric SPK position/velocity records remain sampled.
/// Nonsidereal output differentiates the full selected frame, while base
/// sidereal output retains its sampled angular drift and analytic range rate.
/// Projected sidereal and LEB/SPK center-offset compositions stay sampled.
/// Apparent rates use a dyadic fourth-order derivative inside every actual
/// central stored/SPK read, including composed center-velocity side probes.
/// Tier changes and optional center absence have their own owned domains.
/// Original one-/six-second reads and errors precede this refinement. Base
/// sidereal output retains its original UT-sampled ayanamsha drift; projected
/// prefixes acquire fresh complete support. Opaque sources stay sampled.
///
/// The context's sidereal selection supplies the base offset or the existing
/// fixed-epoch/ECL_T0/SSY_PLANE rigid rotations. Unsupported endpoints,
/// unavailable configuration/data, coverage and numerical failures propagate
/// as typed errors. `skip_closed` controls only closed bound LEB images.
#[allow(clippy::too_many_arguments)]
pub fn calc_pctr_in_context(
    ctx: &EphemerisContext,
    model: crate::time::DeltaTModel,
    iers: bool,
    skip_closed: bool,
    tt: f64,
    target: i32,
    center: i32,
    flags: i32,
) -> Result<(PositionResult, i32)> {
    calc_pctr_in_context_with_anchors(
        ctx,
        model,
        iers,
        skip_closed,
        tt,
        target,
        center,
        flags,
        None,
    )
}

/// PCTR with an optional live anchor owner at its final sidereal reduction.
/// Target/observer, body SPEED samples and actual Time are acquired first.
/// Auto's independent anchor provider is supplied by the caller, never by
/// changing the body provider or reusing its configured model TT.
#[allow(clippy::too_many_arguments)]
pub fn calc_pctr_in_context_with_anchors(
    ctx: &EphemerisContext,
    model: crate::time::DeltaTModel,
    iers: bool,
    skip_closed: bool,
    tt: f64,
    target: i32,
    center: i32,
    flags: i32,
    live_anchors: Option<&dyn crate::sidereal::DynamicAnchorProvider>,
) -> Result<(PositionResult, i32)> {
    if !tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: tt });
    }
    for id in [target, center] {
        if !is_pctr_supported_body(id) {
            return Err(Error::UnknownBody { body_id: id });
        }
    }
    if ctx.mode() == ProviderMode::NativeJpl {
        let provider = super::routing::LocalSpkStateProvider::new(ctx.spk_registry(), ctx.mode())?;
        provider.validate_scalar_tt(target, tt)?;
        let delta_cfg = ctx.delta_t_config(model, iers);
        let reads = RefCell::new(PredictedSpkReads::pctr(false));
        let metadata = |id, hi, lo| provider.pctr_record_offsets_tdb(id, TwoPartEpoch { hi, lo });
        let retarded_tdb = |id, epoch| {
            reads.borrow().query(id, epoch, &metadata, || {
                provider
                    .state_tdb_two_part(id, epoch)
                    .map(PctrState::spk_tdb)
            })
        };
        let state_at = |id, epoch| {
            pctr_track_tt_read(&reads.borrow(), id, epoch, true, &metadata, || {
                provider
                    .pctr_state_two_part(id, epoch)
                    .map(|(state, derivative)| {
                        PctrState::spk_position_tangent((state.pos(), state.vel(), derivative))
                    })
            })
        };
        pctr_context_provider(
            &provider,
            live_anchors,
            ctx,
            &delta_cfg,
            &state_at,
            Some(&retarded_tdb),
            Some(&reads),
            tt,
            target,
            center,
            flags,
        )
        .map_err(|error| provider.public_range_error(target, tt, error))
    } else {
        context_serving::with_sealed_serving_stack(ctx, model, iers, skip_closed, |stack| {
            if let Some((start_jd, end_jd)) = stack.serving.provider().body_coverage(target) {
                if tt < start_jd || tt > end_jd {
                    return Err(Error::EphemerisRange {
                        requested_jd: tt,
                        start_jd,
                        end_jd,
                        body: Some(target),
                    });
                }
            }
            let provider = stack.serving.provider();
            #[cfg(test)]
            let _observed_provider =
                observations::PctrProviderGuard::enter(provider, ctx.leb_files());
            let source = PctrSource::from_leb_provider(provider);
            let reads = RefCell::new(PredictedSpkReads::pctr(false));
            let state_at = |id, epoch: TwoPartEpoch| {
                context_serving::pctr_state_with_owned_support(provider, &reads.borrow(), id, epoch)
                    .map(|state| source.state(state))
                    .map_err(pctr_public_leb_error)
            };
            pctr_context_provider(
                provider,
                live_anchors,
                ctx,
                stack.delta_cfg,
                &state_at,
                None,
                Some(&reads),
                tt,
                target,
                center,
                flags,
            )
        })
    }
}

/// Nonsidereal planet-centric calculation over retained UT1/TT parts.
///
/// `whole` and `tt_fraction` are reception TT days; `original_ut1_fraction`
/// retains the corresponding UT1. Native SPK preserves observation TDB parts;
/// stored LEB evaluates its scalar TT channel. HELCTR/BARYCTR do not change
/// the explicit observer here: TRUEPOS, NOGDEFL and NOABERR select corrections.
/// Final XYZ/RADIANS shaping and implied returned flags are applied once.
///
/// The enclosing controller validates the original epoch and acquires data.
/// This reducer performs no file loading or surface-site reduction. State,
/// frame and coverage errors propagate; a promoted Type-2 coverage error names
/// the original UT1 request, while other provider errors retain their type.
#[allow(clippy::too_many_arguments)]
pub fn calc_pctr_context_ut1_parts_in_context(
    ctx: &EphemerisContext,
    model: crate::time::DeltaTModel,
    iers: bool,
    skip_closed: bool,
    whole: f64,
    tt_fraction: f64,
    original_ut1_fraction: f64,
    target: i32,
    center: i32,
    calc_flags: i32,
    normalized_flags: i32,
) -> Result<(PositionResult, i32)> {
    for body in [target, center] {
        if !is_pctr_supported_body(body) {
            return Err(Error::UnknownBody { body_id: body });
        }
    }
    let epoch = TwoPartEpoch {
        hi: whole,
        lo: tt_fraction,
    };
    let original_ut = whole + original_ut1_fraction;
    let offset = |_: f64, _: i32| {
        Err(Error::ConfigurationMissing {
            missing_config: "context_pctr_sidereal_delegate",
            suggestion: "sidereal Context PCTR must use the original module TT delegation",
        })
    };
    let finish = |state_at: &PctrStateAt<'_>,
                  tdb: Option<&PctrStateAt<'_>>,
                  aliases: &dyn Fn() -> Result<()>,
                  reads: &PredictedSpkReads|
     -> Result<(PositionResult, i32)> {
        let value = pctr_vector_result_tracked(
            state_at,
            tdb,
            &offset,
            epoch,
            target,
            center,
            calc_flags,
            Some(aliases),
            PctrRatePolicy::from_flags(normalized_flags),
            Some(reads),
        )?;
        Ok(value.finalize(calc_flags, normalized_flags))
    };
    if ctx.mode() == ProviderMode::NativeJpl {
        let provider = super::routing::LocalSpkStateProvider::new(ctx.spk_registry(), ctx.mode())?;
        let reads = PredictedSpkReads::pctr(pctr_owned_rates_enabled(normalized_flags));
        let metadata = |id, hi, lo| provider.pctr_record_offsets_tdb(id, TwoPartEpoch { hi, lo });
        let aliases = || {
            provider.validate_source_alias(target)?;
            provider.validate_source_alias(center)
        };
        let state_at = |body, epoch| {
            pctr_track_tt_read(&reads, body, epoch, false, &metadata, || {
                provider
                    .pctr_source_time_state(body, epoch)
                    .map(PctrState::spk_position_tangent)
            })
        };
        let tdb = |body, epoch| {
            reads.query(body, epoch, &metadata, || {
                provider
                    .pctr_source_tdb_state(body, epoch)
                    .map(PctrState::spk_tdb)
            })
        };
        finish(&state_at, Some(&tdb), &aliases, &reads)
            .map_err(|error| provider.pctr_range_error(target, original_ut, error))
    } else {
        context_serving::with_sealed_serving_stack(ctx, model, iers, skip_closed, |stack| {
            let provider = stack.serving.provider();
            #[cfg(test)]
            let _observed_provider =
                observations::PctrProviderGuard::enter(provider, ctx.leb_files());
            let source = PctrSource::from_leb_provider(provider);
            let reads = PredictedSpkReads::pctr(pctr_owned_rates_enabled(normalized_flags));
            let aliases = || {
                for body in [target, center] {
                    if !provider.has_body(body) {
                        return Err(pctr_missing_vector_alias(body));
                    }
                }
                Ok(())
            };
            let state_at = |body, epoch: TwoPartEpoch| {
                context_serving::pctr_state_with_owned_support(provider, &reads, body, epoch)
                    .map(|state| source.state(state))
                    .map_err(|error| provider.vector_error(error))
            };
            finish(&state_at, None, &aliases, &reads)
        })
    }
}

/// A requested body has no vector in the selected LEB provider.
pub(crate) fn pctr_missing_vector_alias(body: i32) -> Error {
    Error::UnknownBody { body_id: body }
}

fn pctr_owned_rates_enabled(flags: i32) -> bool {
    flags & (FLG_SPEED | FLG_TRUEPOS) == FLG_SPEED
}

/// Track actual TDB operands, rather than treating TT offsets as record offsets.
/// Clock errors still come from the original state read in its original order.
fn pctr_track_tt_read<T>(
    reads: &PredictedSpkReads,
    body: i32,
    epoch: TwoPartEpoch,
    checked_clock: bool,
    metadata: &dyn Fn(i32, f64, f64) -> Result<(f64, f64)>,
    read: impl FnOnce() -> Result<T>,
) -> Result<T> {
    reads.before_read(epoch.value())?;
    let correction = if checked_clock {
        match reference_frames::tdb_minus_tt_seconds(epoch.hi, epoch.lo) {
            Ok(value) => value,
            Err(_) => return read(),
        }
    } else {
        crate::time::source_tdb_minus_tt_seconds(epoch.hi, epoch.lo)
    };
    let tdb = TwoPartEpoch {
        hi: epoch.hi,
        lo: epoch.lo + correction / 86400.0,
    };
    reads.query(body, tdb, metadata, read)
}

/// Bind independent vector-state and sidereal clocks to the public PCTR entry.
/// LEB consumes scalar TT at its reader boundary; SPK retains both parts.
/// Explicit TDB overrides are read directly when the provider supplies them.
#[allow(clippy::too_many_arguments)]
fn pctr_context_provider<P: StateProvider>(
    provider: &P,
    live_anchors: Option<&dyn crate::sidereal::DynamicAnchorProvider>,
    ctx: &EphemerisContext,
    delta_cfg: &crate::time::DeltaTConfig<'_>,
    state_at: &PctrStateAt<'_>,
    retarded_tdb: Option<&PctrStateAt<'_>>,
    reads: Option<&RefCell<PredictedSpkReads>>,
    tt: f64,
    target: i32,
    center: i32,
    flags: i32,
) -> Result<(PositionResult, i32)> {
    // Public range validation precedes alias resolution. Both endpoint
    // aliases then resolve before the first state read, including when the
    // other endpoint is present but its derived sample is out of range.
    for body_id in [target, center] {
        if !provider.has_body(body_id) {
            return Err(Error::UnknownBody { body_id });
        }
    }
    let (selection, _) = crate::houses::sidereal_selection_for_context(ctx);
    let delta = |ut| crate::time::deltat_with_config(ut, delta_cfg);
    let anchors = crate::stars::SealedAnchorProvider::new(provider);
    let anchors = if let Some(live) = live_anchors {
        anchors.with_dynamic_override(live)
    } else {
        anchors
    };
    let aya = crate::sidereal::SealedAyanamshaServing::new(selection, &delta, Some(&anchors));
    let offset = |ut, calc_flags| {
        if live_anchors.is_none() {
            return aya.for_flags_ut(ut, calc_flags).map(|(v, _)| v);
        }
        // _get_ayanamsa_for_flags owns a fresh global Time after body rates.
        let tt = crate::time::default_timescale_tt_from_ut1(ut)?;
        let (noaberr, nogdefl, use_true) = crate::sidereal::resolve_ayanamsha_flags(calc_flags);
        if use_true {
            aya.true_tt_with(tt, noaberr, nogdefl)
        } else {
            aya.mean_tt_with(tt, noaberr, nogdefl)
        }
        .map(|(value, _)| value)
    };
    let mut mapped = flags;
    if mapped & (FLG_HELCTR | FLG_BARYCTR) != 0 {
        mapped = (mapped | FLG_NOABERR | FLG_NOGDEFL) & !(FLG_HELCTR | FLG_BARYCTR);
    }
    let raw = |jd, _id, calc_flags| {
        let normal = exclusive_ephemeris_bit(calc_flags, false);
        let core_flags = strip_output_flags(normal);
        // Every projected sidereal prefix is a separate physical place.
        // Its support must come from that call's complete central reads.
        let active_reads = reads.map(|slot| {
            *slot.borrow_mut() = PredictedSpkReads::pctr(pctr_owned_rates_enabled(core_flags));
            slot.borrow()
        });
        // A projection subrequest removes SIDEREAL and later reconstructs
        // the full Cartesian velocity. Preserve that complete sampled route;
        // qualifying its radial component alone would not bound angular rates.
        let rate_policy = if mapped & FLG_SIDEREAL != 0 && core_flags & FLG_SIDEREAL == 0 {
            PctrRatePolicy::Sampled
        } else {
            PctrRatePolicy::from_flags(flags)
        };
        let value = pctr_vector_result_tracked(
            state_at,
            retarded_tdb,
            &offset,
            pctr_tt_parts(jd, 0.0),
            target,
            center,
            core_flags,
            None,
            rate_policy,
            active_reads.as_deref(),
        )?;
        Ok(value.finalize(core_flags, normal))
    };
    use crate::sidereal::SiderealRoute;
    let route =
        crate::sidereal::resolve_sidereal_route(mapped, selection.full().0, selection.bits());
    if matches!(route, None | Some(SiderealRoute::BaseCorrection)) {
        return raw(tt, target, mapped);
    }
    // The shared serving owns only the established rigid rotations here.
    // PCTR's UT-sampled ayanamsha drift stays in its own base branch above.
    let mean_tt = |jd, noaberr, nogdefl| aya.mean_tt_with(jd, noaberr, nogdefl).map(|(v, _)| v);
    let mean_ut = |ut| {
        if live_anchors.is_some() {
            aya.mean_tt(crate::time::default_timescale_tt_from_ut1(ut)?)
                .map(|(value, _)| value)
        } else {
            aya.mean_ut(ut).map(|(value, _)| value)
        }
    };
    let dpsi = |jd| {
        let (psi, _) = reference_frames::iau2000a_tenth_microarcsec(jd)?;
        Ok((psi * reference_frames::FRAME_TENTH_USEC_2_RAD).to_degrees())
    };
    let dpsi_rate = |jd| nutation_rate_deg_day(jd, |t| dpsi(t).map(f64::to_radians));
    let (_, t0, _) = selection.full();
    let user_t0 = if selection.mode() == crate::constants::SIDM_USER
        && selection.bits() & crate::constants::SIDBIT_USER_UT != 0
    {
        t0 + delta(t0)?
    } else {
        t0
    };
    let valens_t0 = if selection.mode() == crate::constants::SIDM_VALENS_MOON {
        aya.valens_epoch_tt()?
    } else {
        0.0
    };
    let sidereal = crate::sidereal::SiderealFastServing::new(
        selection, &raw, &mean_tt, &mean_ut, &dpsi, &dpsi_rate, user_t0, valens_t0,
    );
    let (value, ret) = sidereal.compute_tt(tt, target, strip_output_flags(mapped))?;
    Ok(finalize_output_flags(
        value,
        echo_request_bits(ret, mapped),
        mapped,
    ))
}

/// Preserve the actual failing scalar LEB acquisition in public PCTR errors.
/// Missing aliases are unknown bodies; declared range misses retain the
/// sampled TT and that body's stored bounds. SPK and corrupt-data errors
/// pass through for the owning provider's separate public error policy.
fn pctr_public_leb_error(error: Error) -> Error {
    match error {
        Error::LebBodyNotFound { body_id } => Error::UnknownBody { body_id },
        Error::LebOutOfRange {
            body_id,
            jd,
            start,
            end,
        } => Error::EphemerisRange {
            requested_jd: jd,
            start_jd: start,
            end_jd: end,
            body: Some(body_id),
        },
        error => error,
    }
}

/// Derivative provenance is separate from the vector's ICRS axes and origin.
/// In particular, an SPK derivative is per TDB day, while composed LEB
/// planet-center offsets can contain a sampled TT derivative.
#[derive(Clone, Copy, PartialEq, Eq)]
enum PctrDerivative {
    StoredIcrsTt,
    PolynomialIcrsTdb,
    SampledIcrsTt,
    SpkTdb,
    Opaque,
}

/// Actual acquired provider, separate from frame metadata and rounded values.
#[derive(Clone, Copy, PartialEq, Eq)]
enum PctrSource {
    StoredLeb,
    CenterAwareLeb,
}

impl PctrSource {
    fn from_leb_provider(provider: &context_serving::BoundProvider<'_>) -> Self {
        match provider {
            context_serving::BoundProvider::Composite(_)
            | context_serving::BoundProvider::Tiered(_) => Self::StoredLeb,
            context_serving::BoundProvider::WithCenters(_, _) => Self::CenterAwareLeb,
        }
    }

    fn state(self, state: GeometricState) -> PctrState {
        let derivative = match state.frame() {
            StateFrame::IcrsBarycentric => PctrDerivative::StoredIcrsTt,
            StateFrame::IcrsSystemBarycenter if self == Self::StoredLeb => {
                PctrDerivative::StoredIcrsTt
            }
            StateFrame::IcrsSystemBarycenter => PctrDerivative::SampledIcrsTt,
            _ => PctrDerivative::Opaque,
        };
        PctrState {
            position: state.pos(),
            velocity: state.vel(),
            derivative,
            position_tangent_tdb: None,
        }
    }
}

/// Acquired endpoint state, preserving the provider's original velocity.
pub(super) struct PctrState {
    position: [f64; 3],
    velocity: [f64; 3],
    derivative: PctrDerivative,
    position_tangent_tdb: Option<[f64; 3]>,
}

impl PctrState {
    fn spk_tdb((position, velocity): ([f64; 3], [f64; 3])) -> Self {
        Self {
            position,
            velocity,
            derivative: PctrDerivative::SpkTdb,
            position_tangent_tdb: None,
        }
    }

    /// Retain a separately certified position tangent without changing stored V.
    fn spk_position_tangent(
        (position, velocity, tangent): ([f64; 3], [f64; 3], Option<[f64; 3]>),
    ) -> Self {
        Self {
            position,
            velocity,
            derivative: if tangent.is_some() {
                PctrDerivative::PolynomialIcrsTdb
            } else {
                PctrDerivative::SpkTdb
            },
            position_tangent_tdb: tangent,
        }
    }
}

impl From<([f64; 3], [f64; 3])> for PctrState {
    fn from((position, velocity): ([f64; 3], [f64; 3])) -> Self {
        Self {
            position,
            velocity,
            derivative: PctrDerivative::Opaque,
            position_tangent_tdb: None,
        }
    }
}

/// Only certified stored-polynomial geometry supplies a complete TT tangent.
/// Base sidereal longitude retains its independently sampled anchor drift;
/// its range alone can use the already qualified stored tangent.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(super) enum PctrRatePolicy {
    Sampled,
    StoredRadial,
    StoredGeometric,
}

impl PctrRatePolicy {
    fn from_flags(flags: i32) -> Self {
        if flags & (FLG_TRUEPOS | FLG_SPEED) != FLG_TRUEPOS | FLG_SPEED {
            Self::Sampled
        } else if flags & FLG_SIDEREAL == 0 {
            Self::StoredGeometric
        } else if flags & FLG_XYZ == 0 {
            Self::StoredRadial
        } else {
            Self::Sampled
        }
    }
}

/// Barycentric vector-state reader at a retained two-part TT or TDB epoch.
type PctrStateAt<'a> = dyn Fn(i32, TwoPartEpoch) -> Result<PctrState> + 'a;

#[cfg(test)]
fn observe_state(role: &observations::PctrRoleGuard, state: &PctrState) {
    let provenance = match state.derivative {
        PctrDerivative::StoredIcrsTt => "stored_icrs_tt",
        PctrDerivative::PolynomialIcrsTdb => "position_polynomial_icrs_tdb",
        PctrDerivative::SampledIcrsTt => "composed_icrs_tt",
        PctrDerivative::SpkTdb => "spk_physical_velocity_tdb",
        PctrDerivative::Opaque => "opaque",
    };
    let bits = |value: f64| format!("{:016x}", value.to_bits());
    role.event(
        serde_json::json!({"kind":"pctr_state", "position_au_bits":state.position.map(bits),
        "physical_velocity_au_day_bits":state.velocity.map(bits), "provenance":provenance,
        "position_tangent_tdb_bits":state.position_tangent_tdb.map(|value|value.map(bits))}),
    );
}

#[cfg(test)]
fn observe_role(
    name: &str,
    body: i32,
    epoch: TwoPartEpoch,
    clock: &str,
) -> observations::PctrRoleGuard {
    observations::PctrRoleGuard::enter(serde_json::json!({"name":name,"body":body,"clock":clock,
        "hi_bits":format!("{:016x}",epoch.hi.to_bits()),"lo_bits":format!("{:016x}",epoch.lo.to_bits())}))
}

struct PctrSample {
    position: [f64; 3],
    radial_rate: Option<f64>,
    tangent: Option<PctrTangent>,
}

/// Central operands retained until all existing side acquisitions succeed.
struct PctrTangent {
    ray: [f64; 3],
    velocity_tt: [f64; 3],
    output_ray: [f64; 3],
    nutation: (f64, f64),
}

/// Retain Cartesian motion through final shaping instead of reconstructing
/// it from rounded angles and undefined pole longitudes.
#[derive(Debug, PartialEq)]
pub(super) struct PctrResult {
    pub(super) values: PositionResult,
    output_velocity: Option<[f64; 3]>,
}

impl PctrResult {
    fn finalize(self, core_flags: i32, output_flags: i32) -> (PositionResult, i32) {
        let (mut values, returned_flags) =
            finalize_output_flags(self.values, core_flags, output_flags);
        if output_flags & FLG_XYZ != 0 {
            if let Some(velocity) = self.output_velocity {
                values[3..].copy_from_slice(&velocity);
            }
        }
        (values, returned_flags)
    }
}

/// Use integer-day/fraction TT for the central observation. Rate samples
/// retain small signed day offsets as a separate component.
pub(super) fn pctr_tt_parts(tt: f64, offset: f64) -> TwoPartEpoch {
    if offset != 0.0 {
        return TwoPartEpoch { hi: tt, lo: offset };
    }
    let whole = tt.floor();
    TwoPartEpoch {
        hi: whole,
        lo: (tt - whole) + offset,
    }
}

/// Planet-centric apparent spherical position in degrees and AU at split TT.
/// Hold the observer at reception while retarding the target for three
/// iterations, then apply enabled deflection and aberration. An explicit
/// TDB callback bypasses the ordinary TT-to-TDB conversion.
fn pctr_vector_position(
    state_at: &PctrStateAt<'_>,
    retarded_tdb: Option<&PctrStateAt<'_>>,
    epoch: TwoPartEpoch,
    target: i32,
    center: i32,
    flags: i32,
    rate_policy: PctrRatePolicy,
) -> Result<PctrSample> {
    pctr_vector_position_acquired(
        state_at,
        retarded_tdb,
        epoch,
        target,
        center,
        flags,
        false,
        rate_policy,
    )
}

/// Planet-centric vector reduction with an optional independently acquired
/// observer state. The default uses the original observer; the supplied
/// callback can reacquire its actual reception state before target retardation.
#[allow(clippy::too_many_arguments)]
fn pctr_vector_position_acquired(
    state_at: &PctrStateAt<'_>,
    retarded_tdb: Option<&PctrStateAt<'_>>,
    epoch: TwoPartEpoch,
    target: i32,
    center: i32,
    flags: i32,
    reread_observer: bool,
    rate_policy: PctrRatePolicy,
) -> Result<PctrSample> {
    #[cfg(test)]
    let _observed_sample = observations::PctrRoleGuard::enter(serde_json::json!({
        "name":"sample", "tt_hi_bits":format!("{:016x}",epoch.hi.to_bits()),
        "tt_lo_bits":format!("{:016x}",epoch.lo.to_bits()), "target":target,"center":center,
        "flags":flags,"reread_observer":reread_observer,"retarded_tdb":retarded_tdb.is_some()}));
    #[cfg(test)]
    let observed_target = observe_role("reception_target", target, epoch, "TT");
    let target_state = state_at(target, epoch)?;
    #[cfg(test)]
    {
        observe_state(&observed_target, &target_state);
        drop(observed_target);
    }
    #[cfg(test)]
    let observed_observer = observe_role("reception_observer", center, epoch, "TT");
    let observer_state = state_at(center, epoch)?;
    #[cfg(test)]
    {
        observe_state(&observed_observer, &observer_state);
        drop(observed_observer);
    }
    let (target_pos, target_vel) = (target_state.position, target_state.velocity);
    let (observer_pos, observer_vel) = (observer_state.position, observer_state.velocity);
    for vector in [target_pos, target_vel, observer_pos, observer_vel] {
        check_finite_vec3(vector)?;
    }
    let mut p = [
        target_pos[0] - observer_pos[0],
        target_pos[1] - observer_pos[1],
        target_pos[2] - observer_pos[2],
    ];
    // Stored LEB channels carry TT slopes. A certified SPK position
    // polynomial carries a TDB slope, so differentiate the actual split
    // source clock once after relative subtraction. Keep the original
    // stored velocities above for apparent-place corrections.
    // Keep the side acquisitions below: their errors remain observable.
    let velocity = if rate_policy != PctrRatePolicy::Sampled
        && flags & (FLG_TRUEPOS | FLG_SPEED) == FLG_TRUEPOS | FLG_SPEED
    {
        match (target_state.derivative, observer_state.derivative) {
            (PctrDerivative::StoredIcrsTt, PctrDerivative::StoredIcrsTt) => {
                Some(std::array::from_fn(|axis| {
                    target_vel[axis] - observer_vel[axis]
                }))
            }
            (PctrDerivative::PolynomialIcrsTdb, PctrDerivative::PolynomialIcrsTdb) => {
                let rate = crate::time::source_tdb_rate(epoch.hi, epoch.lo);
                target_state
                    .position_tangent_tdb
                    .zip(observer_state.position_tangent_tdb)
                    .map(|(target, observer)| {
                        std::array::from_fn(|axis| (target[axis] - observer[axis]) * rate)
                    })
            }
            _ => None,
        }
    } else {
        None
    };
    let radial_rate = velocity.map(|velocity| {
        crate::positions::planet_range::RangeState::from_cartesian(p, velocity).rate
    });
    if flags & FLG_TRUEPOS == 0 {
        let (observer_pos, observer_vel) = if reread_observer {
            #[cfg(test)]
            let observed = observe_role("observer_reread", center, epoch, "TT");
            let state = state_at(center, epoch)?;
            #[cfg(test)]
            observe_state(&observed, &state);
            (state.position, state.velocity)
        } else {
            (observer_pos, observer_vel)
        };
        let tdb_fraction = crate::time::tt_to_tdb_fraction(epoch.hi, epoch.lo)?;
        #[cfg(test)]
        observations::pctr_event(serde_json::json!({"kind":"pctr_clock",
            "tt_hi_bits":format!("{:016x}",epoch.hi.to_bits()),
            "tt_lo_bits":format!("{:016x}",epoch.lo.to_bits()),
            "tdb_lo_bits":format!("{:016x}",tdb_fraction.to_bits())}));
        let mut light_time = 0.0;
        for _ in 0..3 {
            light_time = vec3_dist(p) / C_AU_PER_DAY_PCTR;
            let fraction = tdb_fraction - light_time;
            #[cfg(test)]
            observations::pctr_event(serde_json::json!({"kind":"pctr_delay",
                "light_time_days_bits":format!("{:016x}",light_time.to_bits()),
                "emission_tdb_lo_bits":format!("{:016x}",fraction.to_bits()),
                "incoming_ray_au_bits":p.map(|v|format!("{:016x}",v.to_bits()))}));
            let retarded = if let Some(at_tdb) = retarded_tdb {
                #[cfg(test)]
                let observed = observe_role(
                    "emission_target",
                    target,
                    TwoPartEpoch {
                        hi: epoch.hi,
                        lo: fraction,
                    },
                    "TDB",
                );
                let state = at_tdb(
                    target,
                    TwoPartEpoch {
                        hi: epoch.hi,
                        lo: fraction,
                    },
                )?;
                #[cfg(test)]
                observe_state(&observed, &state);
                state
            } else {
                // Convert the emitted TDB epoch to TT before a scalar stored-state lookup.
                let tt_fraction =
                    fraction - crate::time::tdb_minus_tt_seconds(epoch.hi, fraction)? / 86400.0;
                #[cfg(test)]
                let observed = observe_role(
                    "emission_target",
                    target,
                    TwoPartEpoch {
                        hi: epoch.hi,
                        lo: tt_fraction,
                    },
                    "TT",
                );
                let state = state_at(
                    target,
                    TwoPartEpoch {
                        hi: epoch.hi,
                        lo: tt_fraction,
                    },
                )?;
                #[cfg(test)]
                observe_state(&observed, &state);
                state
            }
            .position;
            p = [
                retarded[0] - observer_pos[0],
                retarded[1] - observer_pos[1],
                retarded[2] - observer_pos[2],
            ];
        }
        if flags & FLG_NOGDEFL == 0 {
            #[cfg(test)]
            observations::pctr_event(serde_json::json!({"kind":"pctr_deflection_input",
                "ray_au_bits":p.map(|v|format!("{:016x}",v.to_bits())),
                "observer_au_bits":observer_pos.map(|v|format!("{:016x}",v.to_bits())),
                "tt_bits":format!("{:016x}",epoch.value().to_bits()),
                "light_time_days_bits":format!("{:016x}",light_time.to_bits())}));
            p = apply_gravitational_deflection_with_lookup(
                p,
                observer_pos,
                epoch.value(),
                light_time,
                |id, jd| {
                    #[cfg(test)]
                    let observed = observe_role("deflector", id, pctr_tt_parts(jd, 0.0), "TT");
                    match state_at(id, pctr_tt_parts(jd, 0.0)) {
                        Ok(state) => {
                            #[cfg(test)]
                            observe_state(&observed, &state);
                            Ok(Some(state.position))
                        }
                        // Omit a missing optional native deflector alias. Actual coverage and
                        // other typed provider failures propagate.
                        Err(Error::SpkUnregisteredBody { .. }) => {
                            #[cfg(test)]
                            observed
                                .event(serde_json::json!({"kind":"pctr_optional_alias_absent"}));
                            Ok(None)
                        }
                        Err(error) => Err(error),
                    }
                },
            )?;
            #[cfg(test)]
            observations::pctr_event(serde_json::json!({"kind":"pctr_deflection_output",
                "ray_au_bits":p.map(|v|format!("{:016x}",v.to_bits()))}));
        }
        if flags & FLG_NOABERR == 0 {
            #[cfg(test)]
            observations::pctr_event(serde_json::json!({"kind":"pctr_aberration_input",
                "ray_au_bits":p.map(|v|format!("{:016x}",v.to_bits())),
                "physical_velocity_au_day_bits":observer_vel.map(|v|format!("{:016x}",v.to_bits())),
                "light_time_days_bits":format!("{:016x}",light_time.to_bits())}));
            p = apply_aberration(p, observer_vel, light_time)?;
            #[cfg(test)]
            observations::pctr_event(serde_json::json!({"kind":"pctr_aberration_output",
                "ray_au_bits":p.map(|v|format!("{:016x}",v.to_bits()))}));
        }
    }
    let tt = epoch.value();
    // PCTR follows the full vector Time nutation carrier, independently of
    // any LEB stored series used by the ordinary fast-calc frame.
    let (psi, deps) = reference_frames::iau2000a_tenth_microarcsec(tt)?;
    let scale = reference_frames::FRAME_TENTH_USEC_2_RAD;
    let frames = frame_matrices_for_date(tt, psi * scale, deps * scale)?;
    let frame_flags = if flags & FLG_J2000 != 0 {
        flags & !FLG_ICRS
    } else {
        flags
    };
    let output_ray = if flags & FLG_NONUT != 0 && flags & (FLG_EQUATORIAL | FLG_J2000) == 0 {
        let matrix = vondrak_precession_matrix(tt, flags & FLG_ICRS == 0)?;
        rotate_equatorial_to_ecliptic(mat3_vec3(matrix, p)?, vondrak_mean_obliquity_rad(tt)?)?
    } else {
        // J2000 output selects the bias-including mean frame before other ICRS flags.
        frame_transform(p, tt, frame_flags, true, frames)?
    };
    let position = cartesian_to_spherical_deg(output_ray)?;
    #[cfg(test)]
    observations::pctr_event(serde_json::json!({"kind":"pctr_frame",
        "tt_bits":format!("{:016x}",tt.to_bits()), "flags":flags,"frame_flags":frame_flags,
        "nutation_radians_bits":[format!("{:016x}",(psi*scale).to_bits()),format!("{:016x}",(deps*scale).to_bits())],
        "input_ray_au_bits":p.map(|v|format!("{:016x}",v.to_bits())),
        "output_ray_au_bits":output_ray.map(|v|format!("{:016x}",v.to_bits())),
        "spherical_degrees_au_bits":position.map(|v|format!("{:016x}",v.to_bits()))}));
    let tangent = if rate_policy == PctrRatePolicy::StoredGeometric {
        velocity.map(|velocity_tt| PctrTangent {
            ray: p,
            velocity_tt,
            output_ray,
            nutation: (psi * scale, deps * scale),
        })
    } else {
        None
    };
    Ok(PctrSample {
        position,
        radial_rate,
        tangent,
    })
}

/// PCTR SPEED and base sidereal tail over split TT vector samples.
/// UT-sampled ayanamsha drift uses the observation epoch's independent UT1.
#[cfg(test)]
#[allow(clippy::too_many_arguments)]
pub(super) fn pctr_vector_result(
    state_at: &PctrStateAt<'_>,
    retarded_tdb: Option<&PctrStateAt<'_>>,
    offset: &dyn Fn(f64, i32) -> Result<f64>,
    tt: f64,
    target: i32,
    center: i32,
    flags: i32,
    rate_policy: PctrRatePolicy,
) -> Result<PctrResult> {
    pctr_vector_result_acquired(
        state_at,
        retarded_tdb,
        offset,
        pctr_tt_parts(tt, 0.0),
        target,
        center,
        flags,
        None,
        rate_policy,
    )
}

/// Shared planet-centric position/rate reduction at retained TT parts.
/// An optional alias callback validates each recursively sampled body before
/// its state is consumed. SPEED evaluates symmetric one-/six-second TT offsets.
#[cfg(test)]
#[allow(clippy::too_many_arguments)]
pub(super) fn pctr_vector_result_acquired(
    state_at: &PctrStateAt<'_>,
    retarded_tdb: Option<&PctrStateAt<'_>>,
    offset: &dyn Fn(f64, i32) -> Result<f64>,
    epoch: TwoPartEpoch,
    target: i32,
    center: i32,
    flags: i32,
    aliases: Option<&dyn Fn() -> Result<()>>,
    rate_policy: PctrRatePolicy,
) -> Result<PctrResult> {
    pctr_vector_result_tracked(
        state_at,
        retarded_tdb,
        offset,
        epoch,
        target,
        center,
        flags,
        aliases,
        rate_policy,
        None,
    )
}

/// Supplementary derivatives may read only the complete original central roles.
#[allow(clippy::too_many_arguments)]
fn pctr_vector_result_tracked(
    state_at: &PctrStateAt<'_>,
    retarded_tdb: Option<&PctrStateAt<'_>>,
    offset: &dyn Fn(f64, i32) -> Result<f64>,
    epoch: TwoPartEpoch,
    target: i32,
    center: i32,
    flags: i32,
    aliases: Option<&dyn Fn() -> Result<()>>,
    rate_policy: PctrRatePolicy,
    reads: Option<&PredictedSpkReads>,
) -> Result<PctrResult> {
    let tt = epoch.value();
    let sample = |epoch, flags| {
        if let Some(admit) = aliases {
            admit()?;
        }
        if aliases.is_some() {
            pctr_vector_position_acquired(
                state_at,
                retarded_tdb,
                epoch,
                target,
                center,
                flags,
                true,
                rate_policy,
            )
        } else {
            pctr_vector_position(
                state_at,
                retarded_tdb,
                epoch,
                target,
                center,
                flags,
                rate_policy,
            )
        }
    };
    let p = sample(epoch, flags)?;
    if let Some(reads) = reads {
        reads.finish_central();
    }
    let mut out = [p.position[0], p.position[1], p.position[2], 0.0, 0.0, 0.0];
    let mut output_velocity = None;
    let half = pctr_speed_half_step_days(target);
    if flags & FLG_SPEED != 0 {
        let mut samples = flags & !(FLG_SPEED | FLG_SPEED3);
        if flags & FLG_EQUATORIAL == 0 {
            samples &= !FLG_SIDEREAL;
        }
        let before = sample(pctr_tt_parts(tt, -half), samples)?.position;
        let after = sample(pctr_tt_parts(tt, half), samples)?.position;
        out[3] = (after[0] - before[0]) / (2.0 * half);
        let threshold = 180.0 / (2.0 * half);
        if out[3] > threshold {
            out[3] -= 360.0 / (2.0 * half);
        } else if out[3] < -threshold {
            out[3] += 360.0 / (2.0 * half);
        }
        out[4] = (after[1] - before[1]) / (2.0 * half);
        if let Some(tangent) = p.tangent {
            // Preserve all existing alias/state/position refusals before the
            // additional pure frame-tangent validation.
            let scale = reference_frames::FRAME_TENTH_USEC_2_RAD;
            let rate = if flags & (FLG_J2000 | FLG_NONUT) == 0 {
                reference_frames::iau2000a_rate_tenth_microarcsec_per_day(tt)?
            } else {
                (0.0, 0.0)
            };
            let (matrix, matrix_rate) =
                super::super::frame_velocity::pctr_output_frame_matrix_rate(
                    tt,
                    flags,
                    tangent.nutation,
                    (rate.0 * scale, rate.1 * scale),
                )?;
            let motion = mat3_vec3(matrix, tangent.velocity_tt)?;
            let rotation = mat3_vec3(matrix_rate, tangent.ray)?;
            let velocity = std::array::from_fn(|axis| motion[axis] + rotation[axis]);
            let angular = super::super::derivative::spherical_angular_velocity_direct(
                tangent.output_ray,
                velocity,
            )?;
            out[3] = angular[0];
            out[4] = angular[1];
            output_velocity = Some(velocity);
        }
        out[5] = p
            .radial_rate
            .unwrap_or_else(|| (after[2] - before[2]) / (2.0 * half));
        if let Some(reads) = reads {
            // All original side-place errors above remain observable first.
            // Each new read is admitted against the corresponding central
            // physical record before I/O, including finite emission clocks.
            if let Some(rates) = reads.rates_with_owned_leb(
                &|at| {
                    sample(
                        TwoPartEpoch {
                            hi: epoch.hi,
                            lo: epoch.lo + (at - tt),
                        },
                        samples,
                    )
                    .map(|value| value.position)
                },
                tt,
                1.0 / 8.0,
                false,
            )? {
                out[3..].copy_from_slice(&rates);
            }
        }
    }
    if flags & FLG_SIDEREAL != 0 && flags & FLG_EQUATORIAL == 0 {
        let ut = crate::time::default_timescale_ut1_from_tt(tt)?;
        out[0] = (out[0] - offset(ut, flags)?).rem_euclid(360.0);
        if flags & FLG_SPEED != 0 {
            let before = offset(ut - half, flags)?;
            let after = offset(ut + half, flags)?;
            let drift = (after - before + 180.0).rem_euclid(360.0) - 180.0;
            out[3] -= drift / (2.0 * half);
        }
    }
    Ok(PctrResult {
        values: out,
        output_velocity,
    })
}

#[cfg(test)]
mod owned_apparent_rate_tests {
    use super::*;
    use std::cell::RefCell;

    const EPOCH: f64 = 2451545.0;
    const VELOCITY: f64 = 0.001;
    const FLAGS: i32 = FLG_SPEED | FLG_J2000 | FLG_EQUATORIAL | FLG_NOGDEFL | FLG_NOABERR;

    fn expected_latitude_rate() -> f64 {
        let matrix = crate::positions::FRAME_BIAS_J2000;
        let [x, y, z] = mat3_vec3(matrix, [1.0, 0.0, 0.0]).unwrap();
        let [vx, vy, vz] = mat3_vec3(matrix, [0.0, 0.0, -VELOCITY]).unwrap();
        let rho2 = x * x + y * y;
        ((rho2 * vz - z * (x * vx + y * vy)) / ((rho2 + z * z) * rho2.sqrt())).to_degrees()
    }

    fn closed_place(seam: bool, fail_legacy: bool) -> (Result<PctrResult>, Vec<(i32, f64)>) {
        let reads = PredictedSpkReads::pctr(true);
        let actual = RefCell::new(Vec::new());
        let metadata = |body, _, _| {
            Ok(if seam && body == 14 {
                (-1e-9, 1.0)
            } else {
                (-1.0, 1.0)
            })
        };
        let state = |body, epoch: TwoPartEpoch| {
            reads.query(body, epoch, &metadata, || {
                let at = epoch.value();
                actual.borrow_mut().push((body, at - EPOCH));
                if fail_legacy && body == 14 && at > EPOCH && at - EPOCH < 0.001 {
                    return Err(Error::MathDomain { value: 73.0 });
                }
                let rate = if seam && at < EPOCH - 1e-9 {
                    VELOCITY * 10.0
                } else {
                    VELOCITY
                };
                Ok(PctrState {
                    position: if body == 0 {
                        [1.0, 0.0, 0.0]
                    } else {
                        [0.0, 0.0, rate * (at - EPOCH)]
                    },
                    // An independent advertised channel is deliberately
                    // unrelated to dP; NOABERR must not turn it into a slope.
                    velocity: [99.0, 88.0, 77.0],
                    derivative: PctrDerivative::SpkTdb,
                    position_tangent_tdb: None,
                })
            })
        };
        let result = pctr_vector_result_tracked(
            &state,
            None,
            &|_, _| Ok(0.0),
            pctr_tt_parts(EPOCH, 0.0),
            0,
            14,
            FLAGS,
            None,
            PctrRatePolicy::Sampled,
            Some(&reads),
        );
        (result, actual.into_inner())
    }

    #[test]
    fn owned_apparent_rate_removes_scalar_jd_one_second_bias() {
        let (result, actual) = closed_place(false, false);
        let result = result.unwrap();
        let expected = expected_latitude_rate();
        assert!(
            (result.values[4] - expected).abs() < expected.abs() * 1e-9,
            "actual={} expected={expected}",
            result.values[4]
        );
        // Central, before and after each retain their original five reads
        // before any supplementary dyadic epoch is acquired.
        assert_eq!(
            &actual[..5].iter().map(|x| x.0).collect::<Vec<_>>(),
            &[0, 14, 0, 0, 0]
        );
        assert_eq!(actual[6].0, 14);
        assert_eq!(actual[11].0, 14);
        assert!(actual[6].1 < 0.0 && actual[11].1 > 0.0);
    }

    #[test]
    fn owned_apparent_derivative_does_not_average_across_a_record_seam() {
        let (result, actual) = closed_place(true, false);
        let expected = expected_latitude_rate();
        assert!((result.unwrap().values[4] - expected).abs() < expected.abs() * 1e-9);
        // The legacy left sample is still read. Supplementary observer
        // samples are confined to the actual central record's right side.
        assert!(actual[6].1 < 0.0);
        assert!(
            actual[15..]
                .iter()
                .filter(|x| x.0 == 14)
                .all(|x| x.1 >= 0.0)
        );
    }

    #[test]
    fn original_side_error_precedes_every_supplementary_read() {
        let (result, actual) = closed_place(false, true);
        assert_eq!(result.unwrap_err(), Error::MathDomain { value: 73.0 });
        assert_eq!(actual.len(), 12);
    }

    #[test]
    fn supplementary_samples_retain_the_context_tt_low_word() {
        let reads = PredictedSpkReads::pctr(true);
        let seen = RefCell::new(Vec::new());
        let low = 2f64.powi(-33);
        let epoch = TwoPartEpoch { hi: EPOCH, lo: low };
        assert_eq!(epoch.value(), EPOCH);
        let state = |body, at: TwoPartEpoch| {
            reads.query(body, at, &|_, _, _| Ok((-1.0, 1.0)), || {
                seen.borrow_mut().push(at);
                Ok(PctrState::from((
                    if body == 0 { [1.0, 0.0, 0.0] } else { [0.0; 3] },
                    [0.0; 3],
                )))
            })
        };
        pctr_vector_result_tracked(
            &state,
            None,
            &|_, _| Ok(0.0),
            epoch,
            0,
            14,
            FLAGS,
            None,
            PctrRatePolicy::Sampled,
            Some(&reads),
        )
        .unwrap();
        let seen = seen.borrow();
        assert!(seen.len() > 15);
        for at in &seen[15..] {
            assert_eq!(at.hi, epoch.hi);
            // Dyadic reception offsets preserve this sub-JD residual exactly.
            let displacement = at.lo - low;
            assert_eq!(at.lo - displacement, low);
        }
    }
}

#[cfg(test)]
mod radial_tests {
    use super::*;
    use std::cell::RefCell;

    const FLAGS: i32 = FLG_TRUEPOS | FLG_SPEED | FLG_EQUATORIAL | FLG_J2000;

    #[test]
    fn uncomposed_system_tangent_comes_from_the_actual_provider() {
        let file = crate::leb::LebFile::open(&crate::test_data::core_leb_path()).unwrap();
        let provider = context_serving::provider_from_images(
            std::iter::once((file.path(), &file)),
            crate::context::Tier::Base,
            false,
        )
        .unwrap();
        let source = PctrSource::from_leb_provider(&provider);
        let raw = provider.geometric_state(5, 2451546.0).unwrap();
        assert_eq!(raw.frame(), StateFrame::IcrsSystemBarycenter);
        assert!(!raw.is_planet_center());
        assert!(source.state(raw).derivative == PctrDerivative::StoredIcrsTt);

        // An empty registry currently leaves the same values and false
        // center indication. It cannot certify a center-aware provider's
        // derivative: future/different requests can acquire an SPK offset.
        let center_aware = context_serving::BoundProvider::WithCenters(
            Box::new(provider),
            crate::spk::SpkRegistry::new(),
        );
        let current = center_aware.geometric_state(5, 2451546.0).unwrap();
        assert_eq!(current, raw);
        let source = PctrSource::from_leb_provider(&center_aware);
        assert!(source.state(current).derivative == PctrDerivative::SampledIcrsTt);

        // Center corrections do not apply to direct coord0 channels.
        let direct = center_aware.geometric_state(0, 2451546.0).unwrap();
        assert_eq!(direct.frame(), StateFrame::IcrsBarycentric);
        assert!(source.state(direct).derivative == PctrDerivative::StoredIcrsTt);
        assert!(
            PctrState::spk_tdb((direct.pos(), direct.vel())).derivative == PctrDerivative::SpkTdb
        );

        let tiered = context_serving::BoundProvider::Tiered(
            crate::leb::TieredIndex::new(vec![(
                crate::leb::LebTier::Base,
                crate::leb::CompositeIndex::from_shared_files([&file]).unwrap(),
            )])
            .unwrap(),
        );
        let state = tiered.geometric_state(5, 2451546.0).unwrap();
        assert_eq!((state.pos(), state.vel()), (raw.pos(), raw.vel()));
        assert_eq!(state.tier(), Some(crate::leb::LebTier::Base));
        assert!(
            PctrSource::from_leb_provider(&tiered)
                .state(state)
                .derivative
                == PctrDerivative::StoredIcrsTt
        );
    }

    #[test]
    fn stored_tt_tangent_avoids_scalar_epoch_and_radius_difference_errors() {
        let tt = 2451545.0;
        for rate in [-1.0e-10, 0.0, 1.0e-10] {
            let reads = RefCell::new(Vec::new());
            let state = |body, epoch: TwoPartEpoch| {
                reads.borrow_mut().push((body, epoch));
                let (position, velocity) = if body == 0 {
                    (
                        [3.0 + rate * (epoch.value() - tt), 0.0, 0.0],
                        [rate, 0.0, 0.0],
                    )
                } else {
                    ([2.0, 0.0, 0.0], [0.0; 3])
                };
                Ok(PctrState {
                    position,
                    velocity,
                    derivative: PctrDerivative::StoredIcrsTt,
                    position_tangent_tdb: None,
                })
            };
            let calculate = |policy| {
                pctr_vector_result(
                    &state,
                    None,
                    &|_, _| panic!("nonsidereal"),
                    tt,
                    0,
                    14,
                    FLAGS,
                    policy,
                )
                .unwrap()
                .values
            };
            let sampled = calculate(PctrRatePolicy::Sampled);
            let original_reads = reads.replace(Vec::new());
            let analytic = calculate(PctrRatePolicy::StoredRadial);
            assert_eq!(analytic[..5], sampled[..5]);
            assert_eq!(*reads.borrow(), original_reads);
            assert_eq!(analytic[5], rate);
            if rate != 0.0 {
                assert!((sampled[5] - rate).abs() > rate.abs() * 1e-6);
            }
        }
    }

    #[test]
    fn other_clocks_keep_the_original_sampled_route() {
        for target_clock in [
            PctrDerivative::StoredIcrsTt,
            PctrDerivative::PolynomialIcrsTdb,
            PctrDerivative::SampledIcrsTt,
            PctrDerivative::SpkTdb,
            PctrDerivative::Opaque,
        ] {
            for center_clock in [
                PctrDerivative::StoredIcrsTt,
                PctrDerivative::PolynomialIcrsTdb,
                PctrDerivative::SampledIcrsTt,
                PctrDerivative::SpkTdb,
                PctrDerivative::Opaque,
            ] {
                for flags in [FLAGS, FLAGS | FLG_XYZ] {
                    let state = |body, _: TwoPartEpoch| {
                        Ok(PctrState {
                            position: if body == 0 {
                                [3.0, 0.0, 0.0]
                            } else {
                                [2.0, 0.0, 0.0]
                            },
                            velocity: if body == 0 {
                                [0.125, 0.0, 0.0]
                            } else {
                                [0.0; 3]
                            },
                            position_tangent_tdb: (if body == 0 {
                                target_clock
                            } else {
                                center_clock
                            } == PctrDerivative::PolynomialIcrsTdb)
                                .then_some(if body == 0 {
                                    [0.125, 0.0, 0.0]
                                } else {
                                    [0.0; 3]
                                }),
                            derivative: if body == 0 {
                                target_clock
                            } else {
                                center_clock
                            },
                        })
                    };
                    let result = pctr_vector_result(
                        &state,
                        None,
                        &|_, _| panic!("nonsidereal"),
                        2451545.0,
                        0,
                        14,
                        flags,
                        PctrRatePolicy::from_flags(flags),
                    )
                    .unwrap()
                    .values;
                    let expected = match (target_clock, center_clock) {
                        (PctrDerivative::StoredIcrsTt, PctrDerivative::StoredIcrsTt) => 0.125,
                        (PctrDerivative::PolynomialIcrsTdb, PctrDerivative::PolynomialIcrsTdb) => {
                            0.125 * crate::time::source_tdb_rate(2451545.0, 0.0)
                        }
                        _ => 0.0,
                    };
                    assert_eq!(result[5], expected);
                }
            }
        }
    }

    #[test]
    fn spk_certificate_preserves_raw_velocity_and_non_geometric_routes() {
        let position = [3.0, 4.0, 12.0];
        let velocity = [0.125, -0.25, 0.5];
        let certified = PctrState::spk_position_tangent((position, velocity, Some(velocity)));
        assert_eq!(certified.position, position);
        assert_eq!(certified.velocity, velocity);
        assert!(certified.derivative == PctrDerivative::PolynomialIcrsTdb);
        let other = PctrState::spk_position_tangent((position, velocity, None));
        assert_eq!(other.position, position);
        assert_eq!(other.velocity, velocity);
        assert!(other.derivative == PctrDerivative::SpkTdb);

        for flags in [
            (FLAGS & !FLG_TRUEPOS) | FLG_NOGDEFL,
            (FLAGS & !FLG_SPEED) | FLG_SPEED3,
            FLAGS | FLG_SIDEREAL | FLG_XYZ,
        ] {
            let calls = RefCell::new(Vec::new());
            let calculate = |certificate: bool| {
                let state_at = |body, epoch: TwoPartEpoch| {
                    calls.borrow_mut().push((body, epoch.hi, epoch.lo));
                    Ok(PctrState::spk_position_tangent((
                        if body == 0 { position } else { [0.0; 3] },
                        if body == 0 { velocity } else { [0.0; 3] },
                        certificate.then_some(if body == 0 { velocity } else { [0.0; 3] }),
                    )))
                };
                pctr_vector_result(
                    &state_at,
                    None,
                    &|_, _| Ok(0.0),
                    2451545.0,
                    0,
                    14,
                    flags,
                    PctrRatePolicy::from_flags(flags),
                )
                .unwrap()
            };
            let before = calculate(false);
            let original_calls = calls.replace(Vec::new());
            assert_eq!(calculate(true), before);
            assert_eq!(*calls.borrow(), original_calls);
        }
    }

    #[test]
    fn geometric_spk_rates_use_position_tangent_instead_of_independent_stored_velocity() {
        let reads = RefCell::new(Vec::new());
        let state = |body, epoch: TwoPartEpoch| {
            reads.borrow_mut().push((body, epoch.hi, epoch.lo));
            Ok(PctrState::spk_position_tangent((
                if body == 0 {
                    [3.0, 4.0, 12.0]
                } else {
                    [0.0; 3]
                },
                if body == 0 {
                    [900.0, -800.0, 700.0]
                } else {
                    [0.0; 3]
                },
                Some([0.0; 3]),
            )))
        };
        let flags = FLAGS | FLG_J2000 | FLG_EQUATORIAL;
        let calculate = |policy| {
            pctr_vector_result(
                &state,
                None,
                &|_, _| panic!("nonsidereal"),
                2451545.0,
                0,
                14,
                flags,
                policy,
            )
            .unwrap()
        };
        let sampled = calculate(PctrRatePolicy::Sampled);
        let original_reads = reads.replace(Vec::new());
        let analytic = calculate(PctrRatePolicy::from_flags(flags));
        assert_eq!(analytic.values, sampled.values);
        assert_eq!(analytic.values[3..], [0.0; 3]);
        assert_eq!(*reads.borrow(), original_reads);
    }

    #[test]
    fn analytic_center_does_not_hide_a_side_acquisition_error() {
        let state = |body, epoch: TwoPartEpoch| {
            if epoch.lo < 0.0 {
                return Err(Error::UnknownBody { body_id: body });
            }
            Ok(PctrState {
                position: [1.0; 3],
                velocity: [0.0; 3],
                derivative: PctrDerivative::StoredIcrsTt,
                position_tangent_tdb: None,
            })
        };
        for policy in [
            PctrRatePolicy::StoredRadial,
            PctrRatePolicy::StoredGeometric,
        ] {
            let result = pctr_vector_result(
                &state,
                None,
                &|_, _| panic!("nonsidereal"),
                2451545.0,
                0,
                14,
                FLAGS,
                policy,
            );
            assert_eq!(result, Err(Error::UnknownBody { body_id: 0 }));
        }
    }

    #[test]
    fn coincident_stored_endpoints_keep_the_zero_origin_convention() {
        let state = |_, _: TwoPartEpoch| {
            Ok(PctrState {
                position: [1.0, -2.0, 3.0],
                velocity: [0.01, -0.02, 0.03],
                derivative: PctrDerivative::StoredIcrsTt,
                position_tangent_tdb: None,
            })
        };
        for policy in [
            PctrRatePolicy::StoredRadial,
            PctrRatePolicy::StoredGeometric,
        ] {
            let result = pctr_vector_result(
                &state,
                None,
                &|_, _| panic!("nonsidereal"),
                2451545.0,
                14,
                14,
                FLAGS,
                policy,
            )
            .unwrap()
            .values;
            assert_eq!(result, [0.0; 6]);
        }
    }

    #[test]
    fn cartesian_geometric_motion_uses_the_static_bias_without_spherical_reconstruction() {
        let tt = 2451545.0;
        let tangent = [1e-5, -2e-5, 3e-5];
        let state = |body, epoch: TwoPartEpoch| {
            let (position, velocity) = if body == 0 {
                (
                    std::array::from_fn(|i| [3.0, 2.0, 3.0][i] + tangent[i] * (epoch.value() - tt)),
                    tangent,
                )
            } else {
                ([2.0, 0.0, 0.0], [0.0; 3])
            };
            Ok(PctrState {
                position,
                velocity,
                derivative: PctrDerivative::StoredIcrsTt,
                position_tangent_tdb: None,
            })
        };
        let flags = FLAGS | FLG_XYZ | FLG_ICRS;
        let calculate = |policy| {
            pctr_vector_result(
                &state,
                None,
                &|_, _| panic!("nonsidereal"),
                tt,
                0,
                14,
                flags,
                policy,
            )
            .unwrap()
            .finalize(flags & !FLG_XYZ, flags)
        };
        let sampled = calculate(PctrRatePolicy::Sampled);
        let analytic = calculate(PctrRatePolicy::StoredGeometric);
        assert_eq!(analytic.0[..3], sampled.0[..3]);
        assert_eq!(analytic.1, sampled.1);
        let expected = mat3_vec3(crate::positions::FRAME_BIAS_J2000, tangent).unwrap();
        assert_eq!(analytic.0[3..], expected);
    }

    #[test]
    fn cartesian_shaping_retains_transverse_motion_at_a_spherical_pole() {
        let spherical = [0.0, 90.0, 1.0, 0.0, 0.0, 0.3];
        let original = finalize_output_flags(spherical, FLG_SPEED, FLG_SPEED | FLG_XYZ);
        let result = PctrResult {
            values: spherical,
            output_velocity: Some([0.1, 0.2, 0.3]),
        }
        .finalize(FLG_SPEED, FLG_SPEED | FLG_XYZ);
        assert_eq!(result.0[..3], original.0[..3]);
        assert_eq!(result.0[3..], [0.1, 0.2, 0.3]);
        assert_eq!(result.1, original.1);
    }
}
