//! Global lunar search preserving the separate candidate and default clocks.
use super::*;

/// Find a global lunar eclipse after or before an original UT1 Julian day.
/// mask0 accepts total, partial and penumbral types; centrality is ignored,
/// annular bits are removed and an annular-only request is rejected before
/// acquisition. Returns character bits and ten UT1 Julian days: maximum,
/// reserved0, partial begin/end, total begin/end, penumbral begin/end and two
/// reserved zeros. Missing phases are0. The default-speed Sun/Moon walk is
/// distinct from the caller-flags Moon-only node probe and the geocentric
/// Cartesian maximum/shadow acquisitions supplied by the ephemeris.
///
/// A forward search first probes the preceding full moon within15 days;
/// a boundary in that preceding full-moon seek is ignored, but a failure in
/// its classification propagates. The main walk converts a consumed coverage
/// boundary into exhaustion, and propagates other errors. Maxima within the
/// epoch margin are skipped; lunations advance25 days. The finite horizon is
/// a backstop, never a timeout or no-event approximation. The actual first
/// default position acquisition owns epoch admission;
/// no data, context, observer or clock is changed by this orchestration.
pub fn lun_eclipse_when_source_via(
    eph: &impl LunarGlobalEphemeris,
    jd_start: f64,
    mask: i32,
    backwards: bool,
) -> Result<(i32, [f64; 10])> {
    let mask = normalize_lunar_search_mask(mask)?;
    if !backwards {
        let mut state = |jd| eph.sun_moon_state(jd);
        let previous = match full_moon_source_via(&mut state, jd_start, true) {
            Ok(jd) => Some(jd),
            Err(error) if is_ephemeris_boundary(&error) => None,
            Err(error) => return Err(error),
        };
        if let Some(previous) = previous {
            if jd_start - previous <= ECLIPSE_PRE_PROBE_WINDOW_DAYS {
                if let Some(result) = check_candidate(eph, previous, mask)? {
                    if result.1[0] > jd_start + ECLIPSE_WHEN_EPOCH_MARGIN_DAYS {
                        return Ok(result);
                    }
                }
            }
        }
    }
    let mut jd = jd_start;
    for _ in 0..ECLIPSE_WHEN_MAX_LUNATIONS {
        let mut state = |ut| eph.sun_moon_state(ut);
        let full_moon = if backwards {
            full_moon_source_via(&mut state, jd, true)
        } else {
            full_moon_source_via(&mut state, jd, false)
        };
        let full_moon = match full_moon {
            Ok(value) => value,
            Err(error) if is_ephemeris_boundary(&error) => break,
            Err(error) => return Err(error),
        };
        let result = match check_candidate(eph, full_moon, mask) {
            Ok(value) => value,
            Err(error) if is_ephemeris_boundary(&error) => break,
            Err(error) => return Err(error),
        };
        if let Some((kind, times)) = result {
            let skip = (!backwards && times[0] <= jd_start + ECLIPSE_WHEN_EPOCH_MARGIN_DAYS)
                || (backwards && times[0] >= jd_start - ECLIPSE_WHEN_EPOCH_MARGIN_DAYS);
            if !skip {
                return Ok((kind, times));
            }
        }
        jd = full_moon
            + if backwards {
                -ECLIPSE_WALK_ADVANCE_DAYS
            } else {
                ECLIPSE_WALK_ADVANCE_DAYS
            };
    }
    Err(Error::EclipseSearchExhausted {
        entry: "lun_eclipse_when",
        detail: format!(
            "No matching lunar eclipse found before the ephemeris boundary searching from JD {jd_start}"
        ),
    })
}

/// Next/previous opposition with provider-owned epoch admission. Consume the
/// initial Sun/Moon state before arithmetic, seed with the mean relative speed,
/// then take at most twenty Newton steps with the 0.1-degree/day speed floor
/// and 1e-5-degree convergence. Return the twentieth estimate if unconverged.
/// Errors stop at the consumed provider stage; NaN/infinite epochs are neither
/// replaced nor preempted by an unrelated orchestration guard.
fn full_moon_source_via(
    state: &mut impl FnMut(f64) -> Result<SunMoonState>,
    jd_start: f64,
    backwards: bool,
) -> Result<f64> {
    let first = state(jd_start)?;
    let elongation = solar_eclipses::fold_180_deg(first.moon_lon_deg - first.sun_lon_deg - 180.0);
    let dt = if backwards {
        if elongation > 0.0 {
            -elongation / SYZYGY_MEAN_RELATIVE_SPEED_DEG_DAY
        } else {
            -(360.0 + elongation) / SYZYGY_MEAN_RELATIVE_SPEED_DEG_DAY
        }
    } else if elongation > 0.0 {
        (360.0 - elongation) / SYZYGY_MEAN_RELATIVE_SPEED_DEG_DAY
    } else {
        -elongation / SYZYGY_MEAN_RELATIVE_SPEED_DEG_DAY
    };
    let mut guess = jd_start + dt;
    for _ in 0..SYZYGY_REFINEMENT_ITERATIONS {
        let at = state(guess)?;
        let diff = solar_eclipses::fold_180_deg(at.moon_lon_deg - at.sun_lon_deg - 180.0);
        if diff.abs() < SYZYGY_CONVERGENCE_DEG {
            return Ok(guess);
        }
        let mut relative_speed = at.moon_speed_deg_day - at.sun_speed_deg_day;
        if relative_speed.abs() < SYZYGY_SPEED_FLOOR_DEG_DAY {
            relative_speed = SYZYGY_FALLBACK_SPEED_DEG_DAY;
        }
        guess -= diff / relative_speed;
    }
    Ok(guess)
}

/// Normalize a global lunar type mask without acquiring positions or state.
/// Central/noncentral bits are ignored, annular bits dropped, and an
/// annular-only mask returns InvalidEclipseSearch. mask0 accepts all three
/// lunar types. Returns the normalized bitmask; other unknown bits are removed.
pub fn normalize_lunar_search_mask(mask: i32) -> Result<i32> {
    super::lunar_eclipses::normalize_lunar_filter_mask(mask)
}

/// Classify one full moon: caller-flags Moon longitude, independent mean node,
/// selenocentric depth minimum, canonical shadow type, then ordered contacts.
/// No Sun scalar is acquired solely for the node prefilter. Geometry errors
/// propagate from the actual consumed stage, including the forward preprobe.
fn check_candidate<E: LunarGlobalEphemeris + ?Sized>(
    eph: &E,
    jd_full_moon: f64,
    mask: i32,
) -> Result<Option<(i32, [f64; 10])>> {
    let moon_lon = eph.candidate_moon_longitude_deg(jd_full_moon)?;
    let node_lon = eph.mean_node_lon_deg(jd_full_moon)?;
    let node_dist = moon_node_distance_deg(node_lon, moon_lon)?;
    if node_dist >= LUNAR_NODE_PREFILTER_DEG {
        return Ok(None);
    }
    let mut depth = |jd| {
        let (moon, sun) = eph.geocentric_xyz_au(jd)?;
        selenocentric_depth_deg(moon, sun)
    };
    let jd_max = refine_lunar_eclipse_maximum(&mut depth, jd_full_moon)?;
    let (kind, _) = eph.lun_shadow(jd_max)?;
    if kind == 0 {
        return Ok(None);
    }
    if !((mask & ECL_TOTAL != 0 && kind & ECL_TOTAL != 0)
        || (mask & ECL_PARTIAL != 0 && kind & ECL_PARTIAL != 0)
        || (mask & ECL_PENUMBRAL != 0 && kind & ECL_PENUMBRAL != 0))
    {
        return Ok(None);
    }
    let mut shadow = |jd| eph.lun_shadow(jd);
    Ok(Some((
        kind,
        lun_eclipse_phase_times_via(jd_max, kind, &mut shadow)?,
    )))
}

#[cfg(test)]
mod tests;
