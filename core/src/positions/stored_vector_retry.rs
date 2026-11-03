//! Body-vector evaluation over the retained scalar LEB reader.
//!
//! This is the vector fallback after a completed fast attempt, not a JPL
//! substitution. Geometry, light time and the rotating observer share a split
//! TT clock. Angular derivatives remain inside the actual LEB polynomials
//! consumed by the central observation; each LEB state uses their scalar TT sum.

use super::*;
use context_serving::ServingStack;

/// One successful geometric read and the polynomial that supplied it.
struct StoredReadDomain {
    body: i32,
    epoch: f64,
    interval: Option<(f64, f64)>,
}

#[derive(Default)]
struct StoredReadDomains {
    ordinary: std::cell::RefCell<Vec<StoredReadDomain>>,
    retarded: std::cell::RefCell<Vec<StoredReadDomain>>,
}

impl StoredReadDomains {
    fn offsets(&self) -> (f64, f64) {
        self.ordinary
            .borrow()
            .iter()
            .chain(self.retarded.borrow().iter())
            .fold(
                (f64::NEG_INFINITY, f64::INFINITY),
                |bounds, read| match read.interval {
                    Some((start, end)) => (
                        bounds.0.max(start - read.epoch),
                        bounds.1.min(end - read.epoch),
                    ),
                    None => bounds,
                },
            )
    }
}

/// Sample-local admission to the actual central physical records.
/// Ordinary observer/deflector roles remain distinct from the iterative target
/// reads. An extra converged target iteration may use only the last central
/// target polynomial. A refused read never contributes to a derivative.
struct StoredReadAdmission<'a> {
    domains: &'a StoredReadDomains,
    ordinary: std::cell::Cell<usize>,
    retarded: std::cell::Cell<usize>,
    last_retarded: std::cell::Cell<Option<(i32, f64)>>,
    refused: std::cell::Cell<bool>,
}

impl<'a> StoredReadAdmission<'a> {
    fn new(domains: &'a StoredReadDomains) -> Self {
        Self {
            domains,
            ordinary: std::cell::Cell::new(0),
            retarded: std::cell::Cell::new(0),
            last_retarded: std::cell::Cell::new(None),
            refused: std::cell::Cell::new(false),
        }
    }

    fn complete(&self) -> bool {
        let final_target_owned = match (
            self.domains.retarded.borrow().last(),
            self.last_retarded.get(),
        ) {
            (None, None) => true,
            (Some(read), Some((body, epoch))) => {
                read.body == body
                    && read
                        .interval
                        .is_none_or(|(start, end)| start <= epoch && epoch <= end)
            }
            _ => false,
        };
        !self.refused.get()
            && self.ordinary.get() == self.domains.ordinary.borrow().len()
            && final_target_owned
    }

    fn check(&self, body: i32, epoch: f64, retarded: bool) -> Result<()> {
        let (reads, cursor) = if retarded {
            (&self.domains.retarded, &self.retarded)
        } else {
            (&self.domains.ordinary, &self.ordinary)
        };
        let index = cursor.get();
        let reads = reads.borrow();
        if let Some(read) = reads
            .get(index)
            .or_else(|| retarded.then(|| reads.last()).flatten())
        {
            if read.body == body
                && read
                    .interval
                    .is_none_or(|(start, end)| start <= epoch && epoch <= end)
            {
                cursor.set(index + 1);
                if retarded {
                    self.last_retarded.set(Some((body, epoch)));
                }
                return Ok(());
            }
        }
        self.refused.set(true);
        // Private control flow: the sampler recognizes only this owned-record
        // refusal. Genuine in-domain provider errors retain their own category.
        Err(Error::MathDomain { value: epoch })
    }

    /// A stable unavailable optional support has no successful central role.
    /// Only loss beyond a centrally consumed record marks domain refusal;
    /// in-domain provider failures keep their actual error and existing policy.
    fn failed_read(&self, body: i32, epoch: f64, retarded: bool) {
        let (reads, index) = if retarded {
            (&self.domains.retarded, self.retarded.get())
        } else {
            (&self.domains.ordinary, self.ordinary.get())
        };
        let reads = reads.borrow();
        if let Some(read) = reads
            .get(index)
            .or_else(|| retarded.then(|| reads.last()).flatten())
        {
            if read.body == body
                && read
                    .interval
                    .is_some_and(|(start, end)| epoch < start || epoch > end)
            {
                self.refused.set(true);
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn stored_place(
    provider: &impl StateProvider,
    vector_error: &dyn Fn(Error) -> Error,
    whole: f64,
    fraction: f64,
    body: i32,
    flags: i32,
    site: Option<(f64, f64, f64)>,
    capture: Option<&StoredReadDomains>,
    admission: Option<&StoredReadAdmission<'_>>,
    strict_support: bool,
) -> Result<[f64; 3]> {
    let state = |id, whole, fraction, retarded| {
        if strict_support && !provider.has_body(id) {
            return Err(super::serving::pctr_missing_vector_alias(id));
        }
        let epoch = whole + fraction;
        let value = match provider.geometric_state(id, epoch) {
            Ok(value) => value,
            Err(error) => {
                if let Some(admission) = admission {
                    admission.failed_read(id, epoch, retarded);
                }
                return Err(vector_error(error));
            }
        };
        if let Some(admission) = admission {
            admission.check(id, epoch, retarded)?;
        }
        if let Some(capture) = capture {
            let read = StoredReadDomain {
                body: id,
                epoch,
                interval: provider.smooth_interval_tt(id, epoch)?,
            };
            (if retarded {
                &capture.retarded
            } else {
                &capture.ordinary
            })
            .borrow_mut()
            .push(read);
        }
        Ok((value.pos(), value.vel()))
    };
    super::routing::stored_vector_place_with_read_roles(
        &state,
        whole,
        fraction,
        body,
        flags & !FLG_SPEED,
        site,
        strict_support,
    )
}

/// Derivative of a held LEB observation on its central owned polynomials.
/// The initial support intersects every actually consumed geometric record;
/// full sample-level admission also bounds nonlinear retarded target epochs.
/// The radial stencil, frame selection and zero-distance behavior are retained.
#[allow(clippy::too_many_arguments)]
pub(super) fn stored_vector_position(
    provider: &impl StateProvider,
    vector_error: &dyn Fn(Error) -> Error,
    tt: f64,
    original_ut: Option<f64>,
    body: i32,
    flags: i32,
    site: Option<(f64, f64, f64)>,
    strict_support: bool,
) -> Result<PositionResult> {
    let (whole, fraction) = super::routing::native_position_central_time(tt, original_ut)?;
    let domains = StoredReadDomains::default();
    let central = stored_place(
        provider,
        vector_error,
        whole,
        fraction,
        body,
        flags,
        site,
        (flags & FLG_SPEED != 0).then_some(&domains),
        None,
        strict_support,
    )?;
    let mut output = [central[0], central[1], central[2], 0.0, 0.0, 0.0];
    if flags & FLG_SPEED == 0 || central[2] == 0.0 {
        return Ok(output);
    }
    let mut support = domains.offsets();
    let topocentric = flags & FLG_TOPOCTR != 0 && flags & (FLG_HELCTR | FLG_BARYCTR) == 0;
    let clock_piece = topocentric
        .then(|| crate::smh2016::bundled_deltat_piece(whole + fraction))
        .transpose()?;
    if topocentric {
        let clock = super::routing::local_jpl_site_clock_support(whole, fraction);
        support.0 = support.0.max(clock.0);
        support.1 = support.1.min(clock.1);
    }
    // The derivative retains exact samples; each new offset independently
    // checks all consumed records before publishing its admission result.
    let sample = |offset: f64| -> Result<Option<[f64; 2]>> {
        if offset == 0.0 {
            return Ok(Some([central[0], central[1]]));
        }
        if let Some(piece) = clock_piece {
            if crate::smh2016::bundled_deltat_piece(whole + (fraction + offset))? != piece {
                return Ok(None);
            }
        }
        let admission = StoredReadAdmission::new(&domains);
        let sampled = stored_place(
            provider,
            vector_error,
            whole,
            fraction + offset,
            body,
            flags,
            site,
            None,
            Some(&admission),
            strict_support,
        );
        let value = if admission.refused.get() {
            None
        } else {
            let coordinates = sampled?;
            admission
                .complete()
                .then_some([coordinates[0], coordinates[1]])
        };
        Ok(value)
    };
    [output[3], output[4]] =
        super::derivative::angular_rates_on_owned_domain(&sample, 1.0 / 64.0, support)?;
    let half = if body == MOON {
        6.0 / 86400.0
    } else {
        1.0 / 86400.0
    };
    let before = stored_place(
        provider,
        vector_error,
        tt,
        -half,
        body,
        flags,
        site,
        None,
        None,
        strict_support,
    )?;
    let after = stored_place(
        provider,
        vector_error,
        tt,
        half,
        body,
        flags,
        site,
        None,
        None,
        strict_support,
    )?;
    output[5] = (after[2] - before[2]) / (2.0 * half);
    Ok(output)
}

#[cfg(test)]
#[path = "stored_vector_retry_tests.rs"]
mod tests;

/// Complete ordinary vector retry. The provider owns the actual selected
/// reader throughout central and derivative samples.
/// Configured Delta-T supplies only ayanamsha defining epochs; observation TT
/// and its UT1 partner belong to the independent bundled Time.
#[allow(clippy::too_many_arguments)]
pub(crate) fn calc_stored_vector_retry(
    stack: &ServingStack<'_, '_>,
    ctx: &EphemerisContext,
    jd: f64,
    body: i32,
    flags: i32,
    is_ut: bool,
    live: Option<&dyn crate::sidereal::DynamicAnchorProvider>,
) -> Result<(PositionResult, i32)> {
    let provider = stack.serving.provider();
    let site = ctx
        .observer()
        .map(|site| (site.lon_deg(), site.lat_deg(), site.height_m()));
    let (selection, _) = crate::houses::sidereal_selection_for_context(ctx);
    let delta = |ut| crate::time::deltat_with_config(ut, stack.delta_cfg);
    let anchors = crate::stars::SealedAnchorProvider::new(provider);
    let anchors = if let Some(live) = live {
        anchors.with_dynamic_override(live)
    } else {
        anchors
    };
    let aya = crate::sidereal::SealedAyanamshaServing::new(selection, &delta, Some(&anchors));
    let to_tt = crate::time::default_timescale_tt_from_ut1;
    let request_tt = if is_ut { to_tt(jd)? } else { jd };
    let ut_at = |tt: f64| {
        if is_ut && tt.to_bits() == request_tt.to_bits() {
            Ok(jd)
        } else {
            crate::time::default_timescale_ut1_from_tt(tt)
        }
    };
    let raw = |tt: f64, id: i32, calc_flags: i32| {
        // Selection of target then observer aliases precedes state reads.
        // Default apparent deflector integer lookups occur later, inside the
        // accepted vector factor, and retain their own lookup context.
        for alias in [
            Some(id),
            if calc_flags & FLG_BARYCTR != 0 {
                None
            } else {
                Some(if calc_flags & FLG_HELCTR != 0 {
                    SUN
                } else {
                    EARTH
                })
            },
        ]
        .into_iter()
        .flatten()
        {
            if !provider.has_body(alias) {
                return Err(super::serving::pctr_missing_vector_alias(alias));
            }
        }
        let mut output = stored_vector_position(
            provider,
            &|error| provider.vector_error(error),
            tt,
            is_ut.then_some(jd),
            id,
            calc_flags,
            site,
            true,
        )?;
        if calc_flags & FLG_SIDEREAL != 0 && calc_flags & FLG_EQUATORIAL == 0 && output[2] != 0.0 {
            let ut = ut_at(tt)?;
            let at = |ut| crate::planetary_moons::moon_spk_ayanamsha(ut, calc_flags, &aya);
            let offset = at(ut)?;
            let speed = calc_flags & FLG_SPEED != 0;
            let half = crate::sidereal::SIDEREAL_SPEED_HALF_STEP_DAYS;
            let (before, after) = if speed {
                (at(ut - half)?, at(ut + half)?)
            } else {
                (offset, offset)
            };
            (output[0], output[3]) = crate::sidereal::apply_sidereal_correction(
                output[0], output[3], offset, before, after, speed,
            )?;
        }
        Ok((output, calc_flags))
    };
    // The shared entry owns the final representation/echo and publication,
    // after the complete body/sidereal operation has succeeded.
    let run = |compute_tt: &dyn Fn(f64, i32, i32) -> Result<(PositionResult, i32)>| {
        let engine = CalcEngine {
            ut1_to_tt_jd: &to_tt,
            compute_tt,
            nutation_obliquity_ut: &ecl_nut_result_ut,
            nutation_obliquity_tt: &|tt, flags| ecl_nut_result(tt, flags, &crate::erfa::nut06a),
            jpl_fallback_active: false,
        };
        super::routing::with_serving_trace_source(body, "Skyfield", false, || {
            if is_ut {
                calc_ut(&engine, jd, body, flags)
            } else {
                calc(&engine, jd, body, flags)
            }
        })
    };
    let route = crate::sidereal::resolve_sidereal_route(flags, selection.mode(), selection.bits());
    if matches!(
        route,
        None | Some(crate::sidereal::SiderealRoute::BaseCorrection)
    ) {
        run(&raw)
    } else {
        // The frame engine strips SIDEREAL before its tropical acquisition.
        context_serving::with_fixed_epoch_engine_with_anchors(stack, ctx, &raw, live, |engine| {
            run(engine.compute_tt)
        })
    }
}
