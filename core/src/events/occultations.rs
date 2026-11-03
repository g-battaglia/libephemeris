//! Occultations implementation for the parent module.
//! Public entry points and types are exposed through the parent API.

use super::*;

/// Reject the degenerate Moon-occults-Moon request.
pub fn reject_moon_self_occultation(body_id: Option<i32>) -> Result<()> {
    if body_id == Some(MOON) {
        return Err(Error::InvalidEclipseSearch {
            entry: "lun_occult",
            detail: "lunar occultation of the Moon itself is undefined (body 1).".to_string(),
        });
    }
    Ok(())
}

/// Expand an occultation-type filter into its full bit set.
///
/// An empty filter means "everything possible for this body"; umbral
/// classes imply both central variants; partial implies non-central.
/// Annular classes exist only for the Sun and are stripped from combined
/// filters of other bodies. `body_label` renders the release `{body}`
/// word of the annular-only refusal (the serving binding formats the
/// public id or star name there).
pub fn normalize_occultation_filter(ecltype: i32, is_sun: bool, body_label: &str) -> Result<i32> {
    let mut wanted = ecltype;
    if wanted == (ECL_PARTIAL | ECL_CENTRAL) {
        return Err(Error::InvalidEclipseSearch {
            entry: "lun_occult_when_glob",
            detail: "central partial eclipses do not exist".to_string(),
        });
    }
    if !is_sun {
        let base_classes = wanted & !(ECL_NONCENTRAL | ECL_CENTRAL);
        if base_classes == ECL_ANNULAR || wanted == ECL_ANNULAR_TOTAL {
            return Err(Error::InvalidEclipseSearch {
                entry: "lun_occult_when_glob",
                detail: format!("annular occultations do not exist for object {body_label}"),
            });
        }
        if wanted & (ECL_ANNULAR | ECL_ANNULAR_TOTAL) != 0 {
            wanted &= !(ECL_ANNULAR | ECL_ANNULAR_TOTAL);
        }
    }
    if wanted == 0 {
        wanted = ECL_TOTAL | ECL_PARTIAL | ECL_NONCENTRAL | ECL_CENTRAL;
        if is_sun {
            wanted |= ECL_ANNULAR | ECL_ANNULAR_TOTAL;
        }
    }
    if wanted & (ECL_TOTAL | ECL_ANNULAR | ECL_ANNULAR_TOTAL) != 0 {
        wanted |= ECL_NONCENTRAL | ECL_CENTRAL;
    }
    if wanted & ECL_PARTIAL != 0 {
        wanted |= ECL_NONCENTRAL;
    }
    Ok(wanted)
}

/// Never-occultable-star predicate over the body's geocentric ecliptic
/// latitude in degrees. Pure comparison — NaN compares false,
/// as in the release and the sibling occultation predicates.
pub fn star_never_occulted(body_ecl_lat_deg: f64) -> bool {
    body_ecl_lat_deg.abs() > OCC_STAR_LAT_LIMIT_DEG
}

/// Coarse occultation latitude gate in degrees.
/// Both latitudes must be finite and their gap must not exceed
/// [`OCC_LAT_GATE_DEG`]. Unknown geometry cannot admit an event.
pub fn conjunction_latitude_admits_occultation(body_lat_deg: f64, moon_lat_deg: f64) -> bool {
    body_lat_deg.is_finite()
        && moon_lat_deg.is_finite()
        && (body_lat_deg - moon_lat_deg).abs() <= OCC_LAT_GATE_DEG
}

/// Apparent disc radii in degrees from the topocentric places:
/// `asin(radius_au / dist_au)` for a finite body (`0.0` for a point
/// source, skipping the `asin` call exactly as the release conditional
/// expression does) and `asin(rmoon_au / moon_dist_au)` for the Moon.
/// Pass [`ECL_RMOON_INNER_AU`] as `rmoon_au` for the inner
/// (disappearance/reappearance) contacts. Non-finite words and
/// non-positive distances are refused.
pub fn occultation_radii_deg(
    body_radius_au: f64,
    body_dist_au: f64,
    rmoon_au: f64,
    moon_dist_au: f64,
) -> Result<(f64, f64)> {
    for v in [body_radius_au, body_dist_au, rmoon_au, moon_dist_au] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    if body_dist_au <= 0.0 {
        return Err(Error::NonFiniteArgument {
            value: body_dist_au,
        });
    }
    if moon_dist_au <= 0.0 {
        return Err(Error::NonFiniteArgument {
            value: moon_dist_au,
        });
    }
    let rbody = if body_radius_au == 0.0 {
        0.0
    } else {
        (body_radius_au / body_dist_au).asin().to_degrees()
    };
    let rmoon = (rmoon_au / moon_dist_au).asin().to_degrees();
    Ok((rbody, rmoon))
}

/// In-progress gate for a fixed instant. Non-finite words are refused.
pub fn is_occultation_in_progress(
    min_sep_deg: f64,
    rbody_deg: f64,
    rmoon_deg: f64,
) -> Result<bool> {
    for v in [min_sep_deg, rbody_deg, rmoon_deg] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    Ok(min_sep_deg <= rbody_deg + rmoon_deg)
}

/// Fixed-instant occultation phase from the minimum separation and the
/// apparent radii, all in degrees: deeper than
/// `rbody - rmoon` is annular (reachable only for a body larger than the
/// Moon — never for a point-source star, whose occultations are total
/// once in progress), within `|rbody - rmoon|` is total, anything else
/// overlapping is partial. The caller gates overlap with
/// [`is_occultation_in_progress`] first. Non-finite words are refused.
pub fn classify_occultation_phase(min_sep_deg: f64, rbody_deg: f64, rmoon_deg: f64) -> Result<i32> {
    for v in [min_sep_deg, rbody_deg, rmoon_deg] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    if min_sep_deg < rbody_deg - rmoon_deg {
        Ok(ECL_ANNULAR)
    } else if min_sep_deg < (rbody_deg - rmoon_deg).abs() {
        Ok(ECL_TOTAL)
    } else {
        Ok(ECL_PARTIAL)
    }
}

/// Outer-contact residual in degrees. The contact solve
/// roots this over the separation provider. Non-finite words are
/// refused.
pub fn occultation_outer_residual_deg(sep_deg: f64, rbody_deg: f64, rmoon_deg: f64) -> Result<f64> {
    for v in [sep_deg, rbody_deg, rmoon_deg] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    Ok(sep_deg - (rbody_deg + rmoon_deg))
}

/// Inner-contact residual in degrees. The contact
/// solve roots this over the separation provider. Non-finite
/// words are refused.
pub fn occultation_inner_residual_deg(sep_deg: f64, rbody_deg: f64, rmoon_deg: f64) -> Result<f64> {
    for v in [sep_deg, rbody_deg, rmoon_deg] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    Ok((rbody_deg - rmoon_deg).abs() - sep_deg)
}

impl OccultWhereOutcome {
    /// Occultation character bits (see the struct docs).
    pub fn retflag(&self) -> i32 {
        self.retflag
    }
    /// 10-slot center-point tuple (see the struct docs).
    pub fn geopos(&self) -> [f64; 10] {
        self.geopos
    }
    /// 20-slot local-circumstance array with `attr[3]` = core-shadow
    /// diameter at the center point, km.
    pub fn attr(&self) -> [f64; 20] {
        self.attr
    }
}

/// Assemble the fixed-instant `lun_occult_where` outcome from its two
/// cores:
/// the character bits and center lon/lat come from the where core, the
/// 20 circumstance slots from the how core evaluated at the center point
/// with the WHERE convention, and `attr[3]` is overwritten with the
/// where-core `dcore[0]`.
pub fn compose_occult_where_outcome(
    where_geom: &WhereGeometry,
    mut how_attr: [f64; 20],
) -> OccultWhereOutcome {
    how_attr[3] = where_geom.dcore()[0];
    let mut geopos = [0.0; 10];
    geopos[0] = where_geom.lon_deg();
    geopos[1] = where_geom.lat_deg();
    OccultWhereOutcome {
        retflag: where_geom.retc(),
        geopos,
        attr: how_attr,
    }
}

/// Global outer-contact residual in km (same operation order as the
/// release `lun_occult_when_glob._f_outer` body: `dc[4] / 2.0 + de_km
/// / dc[5] - dc[2]` over the where-core `dcore` words with
/// [`OCC_EARTH_EQ_RADIUS_KM`]). The contact solve (below) roots this
/// over the where-core provider.
pub fn occult_global_outer_residual_km(eph: &impl OccultGlobalEphemeris, jd: f64) -> Result<f64> {
    let dc = eph.where_core(jd)?.dcore();
    Ok(dc[4] / 2.0 + OCC_EARTH_EQ_RADIUS_KM / dc[5] - dc[2])
}

/// Global totality-contact residual in km (same operation order as
/// the release `lun_occult_when_glob._f_total` body: `abs(dc[3]) / 2.0 +
/// de_km / dc[6] - dc[2]`).
pub fn occult_global_total_residual_km(eph: &impl OccultGlobalEphemeris, jd: f64) -> Result<f64> {
    let dc = eph.where_core(jd)?.dcore();
    Ok(dc[3].abs() / 2.0 + OCC_EARTH_EQ_RADIUS_KM / dc[6] - dc[2])
}

/// Global center-line-contact residual in km (same operation order
/// as the release `lun_occult_when_glob._f_central` body: `de_km / dc[6] -
/// dc[2]`).
pub fn occult_global_central_residual_km(eph: &impl OccultGlobalEphemeris, jd: f64) -> Result<f64> {
    let dc = eph.where_core(jd)?.dcore();
    Ok(OCC_EARTH_EQ_RADIUS_KM / dc[6] - dc[2])
}

/// Next global lunar occultation from `jd_start` (same orchestration
/// order as the release `lun_occult_when_glob`, over
/// [`OccultGlobalEphemeris`]): conjunction-by-conjunction walk with the
/// geocentric latitude gate, golden-section maximum refinement, the
/// where-core gate, the epoch-margin direction invariant, the
/// unwanted-type arms, outer/total/central contact solves with the
/// narrow-then-wide fallback, the annular-total sign flip and the
/// right-ascension noon-transit slot.
///
/// `body_id` is the already-resolved numeric id (`None` for a star
/// name); `ecltype` is the raw public filter expanded here via
/// [`normalize_occultation_filter`] (impossible filters are
/// [`Error::InvalidEclipseSearch`]); `is_sun` selects the solar annular
/// classes and `body_label` renders the refusal detail. `backward`
/// selects the search direction and `one_try` the single-conjunction
/// probe (a miss returns `(0, tret)` with the continuation date in
/// `tret[0]`, exactly as in the release).
///
/// Returns `(retflag, tret)` with the 10-slot global phase layout
/// (`[0]` maximum, `[1]` noon-transit, `[2]/[3]` begin/end,
/// `[4]/[5]` totality begin/end, `[6]/[7]` center-line begin/end). A
/// walk that steps off the ephemeris is
/// [`Error::EclipseSearchExhausted`] (never a fabricated event); any
/// other provider failure propagates.
#[allow(clippy::too_many_arguments)]
pub fn lun_occult_when_glob_via(
    eph: &impl OccultGlobalEphemeris,
    jd_start: f64,
    body_id: Option<i32>,
    ecltype: i32,
    is_sun: bool,
    body_label: &str,
    backward: bool,
    one_try: bool,
) -> Result<(i32, [f64; 10])> {
    if !jd_start.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_start });
    }
    lun_occult_when_glob_walk(
        eph, jd_start, body_id, ecltype, is_sun, body_label, backward, one_try,
    )
}
/// Global numeric occultation search over a selected native provider.
/// Reuses the shared conjunction, overlap minimum and contact solves.
#[allow(clippy::too_many_arguments)]
pub(crate) fn lun_occult_when_glob_source_via(
    eph: &impl OccultGlobalEphemeris,
    jd_start: f64,
    body_id: i32,
    ecltype: i32,
    backward: bool,
    one_try: bool,
) -> Result<(i32, [f64; 10])> {
    lun_occult_when_glob_walk(
        eph,
        jd_start,
        Some(body_id),
        ecltype,
        body_id == SUN,
        &body_id.to_string(),
        backward,
        one_try,
    )
}

/// Global occultation search for a named star or resolved numeric target.
/// Named stars use the initial ecliptic-latitude impossibility test.
/// Native data-range errors end the search; other failures propagate.
#[allow(clippy::too_many_arguments)]
pub(crate) fn lun_occult_when_glob_star_source_via(
    eph: &impl OccultGlobalEphemeris,
    jd_start: f64,
    body_id: Option<i32>,
    ecltype: i32,
    body_label: &str,
    backward: bool,
    one_try: bool,
) -> Result<(i32, [f64; 10])> {
    lun_occult_when_glob_walk(
        eph,
        jd_start,
        body_id,
        ecltype,
        body_id == Some(SUN),
        body_label,
        backward,
        one_try,
    )
}

/// Shared bounded conjunction, maximum, classification and contact walk.
#[allow(clippy::too_many_arguments)]
fn lun_occult_when_glob_walk(
    eph: &impl OccultGlobalEphemeris,
    jd_start: f64,
    body_id: Option<i32>,
    ecltype: i32,
    is_sun: bool,
    body_label: &str,
    backward: bool,
    one_try: bool,
) -> Result<(i32, [f64; 10])> {
    reject_moon_self_occultation(body_id)?;
    let ifltype = normalize_occultation_filter(ecltype, is_sun, body_label)?;
    let direction: i32 = if backward { -1 } else { 1 };
    let is_star = body_id.is_none();

    let mut t = jd_start;
    for _ in 0..OCCULT_MAX_CONJUNCTIONS {
        let sought = if is_star {
            super::numeric_occult_global::source_conjunction_with_star(
                &mut |jd| eph.body_lonlat_deg(jd),
                &mut |jd| eph.moon_lonlat_deg(jd),
                t,
                direction,
                Some(body_label),
            )?
        } else {
            super::numeric_occult_global::source_conjunction(
                &mut |jd| eph.body_lonlat_deg(jd),
                &mut |jd| eph.moon_lonlat_deg(jd),
                t,
                direction,
            )?
        };
        let (t_conj, blat, moon_lat) = match sought {
            None => break,
            Some(found) => found,
        };

        // Latitude gate: no occultation possible this conjunction.
        if !conjunction_latitude_admits_occultation(blat, moon_lat) {
            if one_try {
                let mut tret = [0.0; 10];
                tret[0] = t_conj + direction as f64;
                return Ok((0, tret));
            }
            t = t_conj + direction as f64 * OCC_CONJUNCTION_HOP_DAYS;
            continue;
        }

        // Time of maximum: deepest limb overlap (geocentric).
        let mut sep = |jd: f64| eph.sep_minus_radii_deg(jd);
        let tjd = golden_section_min(
            &mut sep,
            t_conj - OCC_MAX_SEARCH_HALF_WINDOW_DAYS,
            t_conj + OCC_MAX_SEARCH_HALF_WINDOW_DAYS,
            GOLDEN_MIN_TOLERANCE_DAYS,
            GOLDEN_SECTION_ITERATIONS,
        )?;

        let where_max = eph.where_core(tjd)?;
        if where_max.retc() == 0 {
            if one_try {
                let mut tret = [0.0; 10];
                tret[0] = tjd;
                return Ok((0, tret));
            }
            t = tjd + direction as f64 * OCC_CONJUNCTION_HOP_DAYS;
            continue;
        }

        let mut tret = [0.0; 10];
        tret[0] = tjd;
        if (backward && tret[0] >= jd_start - ECLIPSE_WHEN_EPOCH_MARGIN_DAYS)
            || (!backward && tret[0] <= jd_start + ECLIPSE_WHEN_EPOCH_MARGIN_DAYS)
        {
            t = tjd + direction as f64 * OCC_CONJUNCTION_HOP_DAYS;
            continue;
        }

        let mut retflag = where_max.retc();
        let unwanted = [ECL_NONCENTRAL, ECL_CENTRAL, ECL_ANNULAR, ECL_PARTIAL]
            .iter()
            .any(|bit| (retflag & bit) != 0 && (ifltype & bit) == 0)
            || ((retflag & ECL_TOTAL) != 0 && (ifltype & (ECL_TOTAL | ECL_ANNULAR_TOTAL)) == 0);
        if unwanted {
            if one_try {
                tret[0] = tjd;
                return Ok((0, tret));
            }
            t = tjd + direction as f64 * OCC_CONJUNCTION_HOP_DAYS;
            continue;
        }

        // Phase times from the shadow geometry: occultation begin/end,
        // totality begin/end, center-line begin/end.
        let mut f_outer = |jd: f64| occult_global_outer_residual_km(eph, jd);
        tret[2] = bisect_sign_change(&mut f_outer, tjd - OCC_GLOBAL_OUTER_WINDOW_DAYS, tjd)?;
        if tret[2] == 0.0 {
            tret[2] = bisect_sign_change(&mut f_outer, tjd - OCC_GLOBAL_OUTER_WIDE_DAYS, tjd)?;
        }
        tret[3] = bisect_sign_change(&mut f_outer, tjd, tjd + OCC_GLOBAL_OUTER_WINDOW_DAYS)?;
        if tret[3] == 0.0 {
            tret[3] = bisect_sign_change(&mut f_outer, tjd, tjd + OCC_GLOBAL_OUTER_WIDE_DAYS)?;
        }
        if (retflag & ECL_PARTIAL) == 0 {
            let mut f_total = |jd: f64| occult_global_total_residual_km(eph, jd);
            tret[4] = bisect_sign_change(&mut f_total, tjd - OCC_GLOBAL_OUTER_WINDOW_DAYS, tjd)?;
            if tret[4] == 0.0 {
                tret[4] = bisect_sign_change(&mut f_total, tjd - OCC_GLOBAL_OUTER_WIDE_DAYS, tjd)?;
            }
            tret[5] = bisect_sign_change(&mut f_total, tjd, tjd + OCC_GLOBAL_OUTER_WINDOW_DAYS)?;
            if tret[5] == 0.0 {
                tret[5] = bisect_sign_change(&mut f_total, tjd, tjd + OCC_GLOBAL_OUTER_WIDE_DAYS)?;
            }
            if (retflag & ECL_NONCENTRAL) == 0 {
                let mut f_central = |jd: f64| occult_global_central_residual_km(eph, jd);
                tret[6] =
                    bisect_sign_change(&mut f_central, tjd - OCC_GLOBAL_OUTER_WINDOW_DAYS, tjd)?;
                if tret[6] == 0.0 {
                    tret[6] =
                        bisect_sign_change(&mut f_central, tjd - OCC_GLOBAL_OUTER_WIDE_DAYS, tjd)?;
                }
                tret[7] =
                    bisect_sign_change(&mut f_central, tjd, tjd + OCC_GLOBAL_OUTER_WINDOW_DAYS)?;
                if tret[7] == 0.0 {
                    tret[7] =
                        bisect_sign_change(&mut f_central, tjd, tjd + OCC_GLOBAL_OUTER_WIDE_DAYS)?;
                }
            }
        }

        // A total occultation that loses its core shadow at begin or end
        // of totality is annular-total (only possible for the Sun).
        if (retflag & ECL_TOTAL) != 0 {
            let d_max = eph.where_core(tret[0])?.dcore()[0];
            let d_beg = if tret[4] != 0.0 {
                eph.where_core(tret[4])?.dcore()[0]
            } else {
                d_max
            };
            let d_end = if tret[5] != 0.0 {
                eph.where_core(tret[5])?.dcore()[0]
            } else {
                d_max
            };
            if d_max * d_beg < 0.0 || d_max * d_end < 0.0 {
                retflag = (retflag & !ECL_TOTAL) | ECL_ANNULAR_TOTAL;
            }
        }

        if ((retflag & ECL_TOTAL) != 0 && (ifltype & ECL_TOTAL) == 0)
            || ((retflag & ECL_ANNULAR_TOTAL) != 0 && (ifltype & ECL_ANNULAR_TOTAL) == 0)
        {
            if one_try {
                tret[0] = tjd;
                return Ok((0, tret));
            }
            t = tjd + direction as f64 * OCC_CONJUNCTION_HOP_DAYS;
            continue;
        }

        // tret[1]: does the occultation happen at local apparent noon
        // somewhere (body-Moon conjunction in right ascension between
        // begin and end)?
        if tret[2] != 0.0 && tret[3] != 0.0 {
            let mut lon_diff = |jd: f64| eph.ra_gap_deg(jd);
            let d_beg = lon_diff(tret[2])?;
            let d_end = lon_diff(tret[3])?;
            tret[1] = if d_beg * d_end >= 0.0 {
                0.0
            } else {
                bisect_sign_change(&mut lon_diff, tret[2], tret[3])?
            };
        }

        return Ok((retflag, tret));
    }

    Err(Error::EclipseSearchExhausted {
        entry: "lun_occult_when_glob",
        detail: format!(
            "no lunar occultation of {body_label} found after JD {jd_start} within the search window"
        ),
    })
}

/// Finish the local occultation attributes (same operation order as
/// local occultation circumstances): the global
/// `ECL_NONCENTRAL` character joins the local `retflag`, `attr[3]`
/// carries the core-shadow width at the global shadow-center point, and
/// the compatibility diameter/disc fractions cap at unity.
///
/// The caps use the `> 1.0` comparison form (not `f64::min`) so a NaN
/// word keeps the release `min(nan, 1.0) == nan` behavior — `f64::min`
/// would return `1.0` there instead.
pub fn finish_occult_local_attr(
    retflag: i32,
    mut how_attr: [f64; 20],
    where_geom: &WhereGeometry,
) -> (i32, [f64; 20]) {
    let retflag = retflag | (where_geom.retc() & ECL_NONCENTRAL);
    how_attr[3] = where_geom.dcore()[0];
    if how_attr[0] > 1.0 {
        how_attr[0] = 1.0;
    }
    if how_attr[2] > 1.0 {
        how_attr[2] = 1.0;
    }
    (retflag, how_attr)
}

/// Next local lunar occultation observable from one place (same
/// local occultation search construction,
/// over [`OccultLocalEphemeris`]): conjunction-by-conjunction walk with the
/// geocentric latitude gate, the topocentric separation-minimum maximum,
/// the in-progress gate, the epoch-margin direction invariant, phase
/// classification, inner contacts (against the smaller lunar radius)
/// and outer contacts with the point-source coincidence for stars,
/// per-contact visibility bits, body rise/set windowing, solar daylight
/// bits and the local-attr assembly.
///
/// `body_id` is the already-resolved numeric id (`None` for a star
/// name); `backward` selects the search direction and `one_try` the
/// single-conjunction probe (a miss returns `(0, tret, zero_attr)` with
/// the continuation date in `tret[0]`, exactly as in the release).
///
/// Returns `(retflag, tret, attr)`: `tret` holds `[0]` maximum,
/// `[1]` first contact, `[2]` second, `[3]` third, `[4]` fourth,
/// `[5]` body rise and `[6]` body set during the occultation (else
/// `0.0`); `attr` holds the local circumstances at the maximum with
/// `attr[3]` = core-shadow width at the global shadow-center point.
/// `retflag` carries the phase bits plus `ECL_VISIBLE` and the
/// per-contact visibility bits, the `ECL_OCC_*_DAYLIGHT` bits and the
/// global `ECL_NONCENTRAL` character. A walk that steps off the
/// ephemeris is [`Error::EclipseSearchExhausted`] (never a fabricated
/// event); any other provider failure propagates.
pub fn lun_occult_when_loc_via(
    eph: &impl OccultLocalEphemeris,
    jd_start: f64,
    body_id: Option<i32>,
    backward: bool,
    one_try: bool,
) -> Result<(i32, [f64; 10], [f64; 20])> {
    if !jd_start.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_start });
    }
    lun_occult_when_loc_walk(eph, jd_start, body_id, backward, one_try, None)
}

/// Local numeric occultation search over a selected native provider.
/// Contacts and visibility use the shared finite-disc search walk.
pub(crate) fn lun_occult_when_loc_source_via(
    eph: &impl OccultLocalEphemeris,
    jd_start: f64,
    body_id: i32,
    backward: bool,
    one_try: bool,
) -> Result<(i32, [f64; 10], [f64; 20])> {
    lun_occult_when_loc_walk(eph, jd_start, Some(body_id), backward, one_try, None)
}

/// Local occultation search for a named point source.
/// The ecliptic-latitude guard precedes the bounded contact/visibility walk.
pub(crate) fn lun_occult_when_loc_star_source_via(
    eph: &impl OccultLocalEphemeris,
    jd_start: f64,
    body_label: &str,
    backward: bool,
    one_try: bool,
) -> Result<(i32, [f64; 10], [f64; 20])> {
    lun_occult_when_loc_walk(eph, jd_start, None, backward, one_try, Some(body_label))
}

/// Shared bounded conjunction, local maximum, contact and visibility walk.
#[allow(clippy::too_many_arguments)]
fn lun_occult_when_loc_walk(
    eph: &impl OccultLocalEphemeris,
    jd_start: f64,
    body_id: Option<i32>,
    backward: bool,
    one_try: bool,
    star_label: Option<&str>,
) -> Result<(i32, [f64; 10], [f64; 20])> {
    reject_moon_self_occultation(body_id)?;
    let is_star = body_id.is_none();
    let direction: i32 = if backward { -1 } else { 1 };
    let zero_attr = [0.0; 20];

    let mut t = jd_start;
    for _ in 0..OCCULT_MAX_CONJUNCTIONS {
        let sought = super::numeric_occult_global::source_conjunction_with_star(
            &mut |jd| eph.body_lonlat_deg(jd),
            &mut |jd| eph.moon_lonlat_deg(jd),
            t,
            direction,
            if is_star {
                Some(star_label.unwrap_or("star"))
            } else {
                None
            },
        )?;
        let (t_conj, blat, moon_lat) = match sought {
            None => break,
            Some(found) => found,
        };

        // Geocentric latitude gate for this conjunction.
        if !conjunction_latitude_admits_occultation(blat, moon_lat) {
            if one_try {
                let mut tret = [0.0; 10];
                tret[0] = t_conj + direction as f64;
                return Ok((0, tret, zero_attr));
            }
            t = t_conj + direction as f64 * OCC_CONJUNCTION_HOP_DAYS;
            continue;
        }

        // Topocentric maximum: deepest limb overlap seen from the place.
        let mut sep_minus = |jd: f64| -> Result<f64> {
            let (rb, rm, sep) = eph.topo_radii_sep_deg(jd, false)?;
            occultation_outer_residual_deg(sep, rb, rm)
        };
        let jd_max = golden_section_min(
            &mut sep_minus,
            t_conj - OCC_MAX_SEARCH_HALF_WINDOW_DAYS,
            t_conj + OCC_MAX_SEARCH_HALF_WINDOW_DAYS,
            GOLDEN_MIN_TOLERANCE_DAYS,
            GOLDEN_SECTION_ITERATIONS,
        )?;

        let (rbody, rmoon, min_sep) = eph.topo_radii_sep_deg(jd_max, false)?;
        if min_sep > rbody + rmoon {
            if one_try {
                let mut tret = [0.0; 10];
                tret[0] = jd_max + direction as f64;
                return Ok((0, tret, zero_attr));
            }
            t = jd_max + direction as f64 * OCC_CONJUNCTION_HOP_DAYS;
            continue;
        }

        if (backward && jd_max >= jd_start - ECLIPSE_WHEN_EPOCH_MARGIN_DAYS)
            || (!backward && jd_max <= jd_start + ECLIPSE_WHEN_EPOCH_MARGIN_DAYS)
        {
            if one_try {
                let mut tret = [0.0; 10];
                tret[0] = jd_max + direction as f64;
                return Ok((0, tret, zero_attr));
            }
            t = jd_max + direction as f64 * OCC_CONJUNCTION_HOP_DAYS;
            continue;
        }

        let phase = classify_occultation_phase(min_sep, rbody, rmoon)?;
        let contact = |f: &mut dyn FnMut(f64) -> Result<f64>, lo, hi| {
            bisect_sign_change(&mut |jd| f(jd), lo, hi)
        };

        // Inner contacts use exact apparent-disc tangency against the
        // smaller inner-contact lunar radius; for point-source stars the
        // outer contacts coincide with them.
        let mut jd_second = 0.0;
        let mut jd_third = 0.0;
        if min_sep < (rbody - rmoon).abs() {
            let mut f_inner = |jd: f64| -> Result<f64> {
                let (rb, rm, sep) = eph.topo_radii_sep_deg(jd, true)?;
                occultation_inner_residual_deg(sep, rb, rm)
            };
            jd_second = contact(&mut f_inner, jd_max - OCC_LOCAL_INNER_WINDOW_DAYS, jd_max)?;
            jd_third = contact(&mut f_inner, jd_max, jd_max + OCC_LOCAL_INNER_WINDOW_DAYS)?;
        }

        let (jd_first, jd_fourth) = if is_star {
            collapse_point_source_contacts(jd_second, jd_third)
        } else {
            let mut f_outer = |jd: f64| -> Result<f64> {
                let (rb, rm, sep) = eph.topo_radii_sep_deg(jd, false)?;
                occultation_outer_residual_deg(sep, rb, rm)
            };
            let mut first = contact(&mut f_outer, jd_max - OCC_LOCAL_OUTER_WINDOW_DAYS, jd_max)?;
            if first == 0.0 {
                first = contact(
                    &mut f_outer,
                    jd_max - 2.0 * OCC_LOCAL_OUTER_WINDOW_DAYS,
                    jd_max,
                )?;
            }
            let mut fourth = contact(&mut f_outer, jd_max, jd_max + OCC_LOCAL_OUTER_WINDOW_DAYS)?;
            if fourth == 0.0 {
                fourth = contact(
                    &mut f_outer,
                    jd_max,
                    jd_max + 2.0 * OCC_LOCAL_OUTER_WINDOW_DAYS,
                )?;
            }
            (first, fourth)
        };

        // Visibility bits from the body's apparent altitude at each
        // contact.
        let mut retflag = phase;
        for (tc, bit) in [
            (jd_max, ECL_MAX_VISIBLE),
            (jd_first, ECL_1ST_VISIBLE),
            (jd_second, ECL_2ND_VISIBLE),
            (jd_third, ECL_3RD_VISIBLE),
            (jd_fourth, ECL_4TH_VISIBLE),
        ] {
            if tc == 0.0 {
                continue;
            }
            if eph.body_alt_deg(tc)? > 0.0 {
                retflag |= ECL_VISIBLE | bit;
            }
        }

        if (retflag & ECL_VISIBLE) == 0 {
            if one_try {
                let mut tret = [0.0; 10];
                tret[0] = jd_max + direction as f64;
                return Ok((0, tret, zero_attr));
            }
            t = jd_max + direction as f64 * OCC_CONJUNCTION_HOP_DAYS;
            continue;
        }

        // Rise and set of the occulted body during the occultation.
        let mut tret = [0.0; 10];
        tret[0] = jd_max;
        tret[1] = jd_first;
        tret[2] = jd_second;
        tret[3] = jd_third;
        tret[4] = jd_fourth;

        let (rise, set) = eph.body_rise_set(jd_first - 0.1)?;
        if let (Some(tr_rise), Some(tr_set)) = (rise, set) {
            if jd_first < tr_rise && tr_rise < jd_fourth {
                tret[5] = tr_rise;
            }
            if jd_first < tr_set && tr_set < jd_fourth {
                tret[6] = tr_set;
            }
        }

        // Daylight flags: is the Sun up at first/fourth contact? (The
        // next sunset coming before the next sunrise means daytime.)
        for (tc, bit) in [
            (jd_first, ECL_OCC_BEG_DAYLIGHT),
            (jd_fourth, ECL_OCC_END_DAYLIGHT),
        ] {
            if tc == 0.0 {
                continue;
            }
            let (sunrise, sunset) = eph.sun_rise_set(tc)?;
            if let (Some(tr_sr), Some(tr_ss)) = (sunrise, sunset) {
                if tr_ss < tr_sr {
                    retflag |= bit;
                }
            }
        }

        // Local circumstances at maximum; global character and core
        // width from the shadow geometry.
        let how_attr = eph.local_how_attr(jd_max)?;
        let where_geom = eph.where_core(jd_max)?;
        let (retflag, attr) = finish_occult_local_attr(retflag, how_attr, &where_geom);

        return Ok((retflag, tret, attr));
    }

    Err(Error::EclipseSearchExhausted {
        entry: "lun_occult_when_loc",
        detail: format!(
            "no observable lunar occultation found within the search window from JD {jd_start}"
        ),
    })
}

impl OccultBodyLeg {
    /// True for the star leg (point source, geocentric throughout).
    pub fn is_star(&self) -> bool {
        matches!(self, OccultBodyLeg::Star { .. })
    }

    /// Numeric id for the radius lookup (`None` = point-source star,
    /// release `_occ_body_radius_au` string arm).
    pub fn body_id_opt(&self) -> Option<i32> {
        match self {
            OccultBodyLeg::Planet(id) => Some(*id),
            OccultBodyLeg::Star { .. } => None,
        }
    }

    /// True when the leg is the Sun (selects the solar annular
    /// classes in the type filter).
    pub fn is_sun(&self) -> bool {
        matches!(self, OccultBodyLeg::Planet(id) if *id == crate::constants::SUN)
    }

    /// Label rendering the refusal detail.
    pub fn label(&self) -> String {
        match self {
            OccultBodyLeg::Planet(id) => id.to_string(),
            OccultBodyLeg::Star { name, .. } => name.clone(),
        }
    }

    /// Disc radius of the occulted body in AU.
    pub fn radius_au(&self) -> f64 {
        occ_body_radius_au(self.body_id_opt())
    }
}

impl<'a, P> SealedOccultGlobalServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    /// Checked constructor for the planet leg: borrows the sealed
    /// serving pair and the Delta-T provider, resolving the
    /// sidereal-continuity offsets from the provider's own boundary
    /// words. Flag selections outside the verified tropical basis are
    /// refused (same gate as the eclipse servings). A negative body id
    /// folds to the Sun. No I/O, no globals.
    pub fn try_new_planet(
        serving: &'a crate::positions::SealedServing<'a, P>,
        star_serving: &'a crate::stars::SealedStarServing<'a, P>,
        deltat_at: &'a dyn Fn(f64) -> Result<f64>,
        calc_flags: i32,
        body_id: i32,
    ) -> Result<Self> {
        let body_id = if body_id < 0 {
            crate::constants::SUN
        } else {
            body_id
        };
        Self::checked(
            serving,
            star_serving,
            deltat_at,
            calc_flags,
            OccultBodyLeg::Planet(body_id),
            "occ_glob_serving",
        )
    }

    /// Checked constructor for the star leg: resolves the star name
    /// through the v1 (`fixstar_ut`) grammar (same resolution as
    /// [`SealedStarEventServing`]).
    /// No I/O, no globals.
    pub fn try_new_star(
        serving: &'a crate::positions::SealedServing<'a, P>,
        star_serving: &'a crate::stars::SealedStarServing<'a, P>,
        deltat_at: &'a dyn Fn(f64) -> Result<f64>,
        calc_flags: i32,
        star_name: &str,
    ) -> Result<Self> {
        let star = match crate::stars::resolve_star_ref(star_name) {
            Ok(matched) => matched,
            Err(crate::stars::StarResolveError::Miss(failure)) => {
                return Err(failure.to_error(star_name));
            }
            Err(crate::stars::StarResolveError::Refused(err)) => return Err(err),
        };
        Self::checked(
            serving,
            star_serving,
            deltat_at,
            calc_flags,
            OccultBodyLeg::Star {
                star,
                name: star_name.to_string(),
            },
            "occ_glob_serving",
        )
    }

    /// Shared checked construction (flag gate plus Delta-T boundary
    /// words).
    fn checked(
        serving: &'a crate::positions::SealedServing<'a, P>,
        star_serving: &'a crate::stars::SealedStarServing<'a, P>,
        deltat_at: &'a dyn Fn(f64) -> Result<f64>,
        calc_flags: i32,
        body: OccultBodyLeg,
        entry: &'static str,
    ) -> Result<Self> {
        use crate::constants::{
            FLG_BARYCTR, FLG_HELCTR, FLG_J2000, FLG_NONUT, FLG_RADIANS, FLG_SIDEREAL, FLG_XYZ,
        };
        if (calc_flags & (FLG_SIDEREAL | FLG_J2000 | FLG_NONUT | FLG_XYZ | FLG_RADIANS)) != 0
            || (calc_flags & (FLG_HELCTR | FLG_BARYCTR)) != 0
        {
            return Err(Error::UnsupportedFlagCombination {
                entry,
                flags: calc_flags,
            });
        }
        let dt_lo = deltat_at(crate::time::SIDEREAL_MODERN_START_JD)?;
        if !dt_lo.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt_lo });
        }
        let dt_hi = deltat_at(crate::time::SIDEREAL_MODERN_END_JD)?;
        if !dt_hi.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt_hi });
        }
        let (offset_lo_deg, offset_hi_deg) =
            crate::time::sidereal_continuity_offsets(dt_lo, dt_hi)?;
        Ok(Self {
            serving,
            star_serving,
            deltat_at,
            calc_flags,
            offset_lo_deg,
            offset_hi_deg,
            body,
        })
    }

    /// Bound body leg.
    pub fn body_leg(&self) -> &OccultBodyLeg {
        &self.body
    }

    /// Bound base calculation flags (without the per-acquisition shaping
    /// bits the release call sites add).
    pub fn calc_flags(&self) -> i32 {
        self.calc_flags
    }

    /// Fixed-instant `lun_occult_where` over the sealed serving (the
    /// occultation ground geometry): the
    /// character bits and shadow-center lon/lat come from the sealed
    /// where core at `jd_ut`, and the 20 circumstance slots from the
    /// sealed how core evaluated at the center point (sea level)
    /// with the WHERE convention, with `attr[3]` overwritten by the
    /// where-core `dcore[0]`. The center site travels explicitly
    /// (explicit-site contract): a temporary center-bound local serving is
    /// built from the same sealed pair, Delta-T source, flags and
    /// body leg, so no global observer state is touched. Numeric
    /// planet leg only (star names stay star-catalog-serving scope); the Moon itself
    /// refuses through the release guard and a negative bound id
    /// already folds to the Sun at construction.
    pub fn lun_occult_where_at_ut(
        &self,
        nutation_at: &dyn Fn(f64) -> Result<(f64, f64)>,
        jd_ut: f64,
    ) -> Result<OccultWhereOutcome> {
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        reject_moon_self_occultation(self.body.body_id_opt())?;
        let where_geom = self.where_geometry_at_ut(jd_ut)?;
        let center = SealedOccultLocalServing::checked(
            self.serving,
            self.star_serving,
            self.deltat_at,
            nutation_at,
            self.calc_flags,
            self.body.clone(),
            where_geom.lon_deg(),
            where_geom.lat_deg(),
            0.0,
            "occ_where_serving",
        )?;
        let how_attr = center.local_how_attr_where(jd_ut)?;
        Ok(compose_occult_where_outcome(&where_geom, how_attr))
    }

    /// TT epoch for a UT epoch through the bound Delta-T provider.
    pub fn jd_tt_for_ut(&self, jd_ut: f64) -> Result<f64> {
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let dt = (self.deltat_at)(jd_ut)?;
        if !dt.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt });
        }
        Ok(jd_ut + dt)
    }

    /// Map a sealed provider miss onto the public range contract
    /// a missing channel inside asset
    /// coverage stays the typed miss, past coverage it names the
    /// requested body with the asset window.
    fn escalate(&self, jd_tt: f64, body_id: i32, err: Error) -> Error {
        match err {
            Error::LebBodyNotFound { .. } | Error::LebOutOfRange { .. } => {
                crate::positions::escalate_sealed_provider_miss(
                    crate::context::ProviderMode::SealedLeb,
                    self.serving.provider(),
                    body_id,
                    jd_tt,
                    err,
                )
            }
            other => other,
        }
    }

    /// Apparent ecliptic words `(lon, lat, dist, dlon, dlat, ddist)` for
    /// one planet id at a TT epoch.
    fn ecliptic_words_at_tt(&self, body_id: i32, jd_tt: f64) -> Result<[f64; 6]> {
        use crate::constants::FLG_SPEED;
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let (pos, _) = crate::positions::serving_compute_tt(
            self.serving,
            jd_tt,
            body_id,
            self.calc_flags | FLG_SPEED,
        )
        .map_err(|e| self.escalate(jd_tt, body_id, e))?;
        Ok(pos)
    }

    /// Geocentric apparent ecliptic `(lon, lat)` of the occulted body,
    /// degrees.
    fn body_lonlat_impl(&self, jd_ut: f64) -> Result<(f64, f64)> {
        match &self.body {
            OccultBodyLeg::Planet(id) => {
                let jd_tt = self.jd_tt_for_ut(jd_ut)?;
                let words = self.ecliptic_words_at_tt(*id, jd_tt)?;
                Ok((words[0], words[1]))
            }
            OccultBodyLeg::Star { star, .. } => {
                use crate::constants::FLG_SPEED;
                let (place, _, _) = self.star_serving.fixstar_ut(star, jd_ut, FLG_SPEED, None)?;
                Ok((place[0], place[1]))
            }
        }
    }

    /// Geocentric apparent ecliptic `(lon, lat)` of the Moon, degrees
    fn moon_lonlat_impl(&self, jd_ut: f64) -> Result<(f64, f64)> {
        use crate::constants::MOON;
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let words = self.ecliptic_words_at_tt(MOON, jd_tt)?;
        Ok((words[0], words[1]))
    }

    /// Geocentric apparent equatorial Cartesian AU vector of the
    /// occulted body.
    /// The planet path runs through the tested
    /// [`crate::positions::apply_output_flags`] spherical-to-Cartesian
    /// conversion of the same sealed RA/Dec/dist words the release
    /// `_cart` helper consumes (same composition as the eclipse
    /// servings).
    fn body_equatorial_xyz_au(&self, jd_ut: f64) -> Result<(f64, f64, f64)> {
        use crate::constants::{FLG_EQUATORIAL, FLG_XYZ};
        match &self.body {
            OccultBodyLeg::Planet(id) => {
                let jd_tt = self.jd_tt_for_ut(jd_ut)?;
                let shape = self.calc_flags | FLG_EQUATORIAL;
                let (pos, _) =
                    crate::positions::serving_compute_tt(self.serving, jd_tt, *id, shape)
                        .map_err(|e| self.escalate(jd_tt, *id, e))?;
                let xyz = crate::positions::apply_output_flags(pos, FLG_XYZ);
                Ok((xyz[0], xyz[1], xyz[2]))
            }
            OccultBodyLeg::Star { star, .. } => {
                let (place, _, _) =
                    self.star_serving
                        .fixstar_ut(star, jd_ut, FLG_EQUATORIAL | FLG_XYZ, None)?;
                Ok((place[0], place[1], place[2]))
            }
        }
    }

    /// Geocentric apparent equatorial Cartesian AU vector of the Moon
    fn moon_equatorial_xyz_au(&self, jd_ut: f64) -> Result<(f64, f64, f64)> {
        use crate::constants::{FLG_EQUATORIAL, FLG_XYZ, MOON};
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let shape = self.calc_flags | FLG_EQUATORIAL;
        let (pos, _) = crate::positions::serving_compute_tt(self.serving, jd_tt, MOON, shape)
            .map_err(|e| self.escalate(jd_tt, MOON, e))?;
        let xyz = crate::positions::apply_output_flags(pos, FLG_XYZ);
        Ok((xyz[0], xyz[1], xyz[2]))
    }

    /// Geocentric separation minus the sum of the apparent radii,
    /// degrees.
    fn sep_minus_radii_impl(&self, jd_ut: f64) -> Result<f64> {
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let body_xyz = self.body_equatorial_xyz_au(jd_ut)?;
        let moon_xyz = self.moon_equatorial_xyz_au(jd_ut)?;
        let db = (body_xyz.0.powi(2) + body_xyz.1.powi(2) + body_xyz.2.powi(2)).sqrt();
        let dm = (moon_xyz.0.powi(2) + moon_xyz.1.powi(2) + moon_xyz.2.powi(2)).sqrt();
        let cosang = (body_xyz.0 * moon_xyz.0 + body_xyz.1 * moon_xyz.1 + body_xyz.2 * moon_xyz.2)
            / (db * dm);
        let sep = cosang.clamp(-1.0, 1.0).acos().to_degrees();
        let rmoon = (ECL_RMOON_AU / dm).asin().to_degrees();
        let body_radius_au = self.body.radius_au();
        let rbody = if body_radius_au != 0.0 {
            (body_radius_au / db).asin().to_degrees()
        } else {
            0.0
        };
        Ok(sep - (rmoon + rbody))
    }

    /// Wrapped body-minus-Moon right-ascension gap, degrees.
    fn ra_gap_impl(&self, jd_ut: f64) -> Result<f64> {
        use crate::constants::{FLG_EQUATORIAL, MOON};
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let body_ra_deg = match &self.body {
            OccultBodyLeg::Planet(id) => {
                let jd_tt = self.jd_tt_for_ut(jd_ut)?;
                let (pos, _) = crate::positions::serving_compute_tt(
                    self.serving,
                    jd_tt,
                    *id,
                    self.calc_flags | FLG_EQUATORIAL,
                )
                .map_err(|e| self.escalate(jd_tt, *id, e))?;
                pos[0]
            }
            OccultBodyLeg::Star { star, .. } => {
                let (place, _, _) =
                    self.star_serving
                        .fixstar_ut(star, jd_ut, FLG_EQUATORIAL, None)?;
                place[0]
            }
        };
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let (moon_pos, _) = crate::positions::serving_compute_tt(
            self.serving,
            jd_tt,
            MOON,
            self.calc_flags | FLG_EQUATORIAL,
        )
        .map_err(|e| self.escalate(jd_tt, MOON, e))?;
        Ok(wrap_ra_gap_deg(body_ra_deg - moon_pos[0]))
    }

    /// Greenwich apparent sidereal time in degrees at a UT epoch
    fn gast_deg_at_ut(&self, jd_ut: f64) -> Result<f64> {
        use crate::time::{SiderealDeltaT, greenwich_apparent_sidereal_time_hours};
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let dt_days = (self.deltat_at)(jd_ut)?;
        if !dt_days.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt_days });
        }
        let jd_tt = jd_ut + dt_days;
        let (dpsi_rad, deps_rad) = (self.serving.nutation_provider())(jd_tt)?;
        let dpsi_deg = dpsi_rad.to_degrees();
        let eps_mean_deg = crate::positions::vondrak_mean_obliquity_deg(jd_tt)?;
        let eps_true_deg = eps_mean_deg + deps_rad.to_degrees();
        let sdt = SiderealDeltaT::try_new(dt_days, self.offset_lo_deg, self.offset_hi_deg)?;
        Ok(greenwich_apparent_sidereal_time_hours(jd_ut, eps_true_deg, dpsi_deg, &sdt)? * 15.0)
    }

    /// Shadow-axis ground geometry at a UT epoch.
    pub(super) fn where_geometry_at_ut(&self, jd_ut: f64) -> Result<WhereGeometry> {
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let moon_xyz = self.moon_equatorial_xyz_au(jd_ut)?;
        let body_xyz = self.body_equatorial_xyz_au(jd_ut)?;
        let gast_deg = self.gast_deg_at_ut(jd_ut)?;
        eclipse_where_core_via(moon_xyz, body_xyz, gast_deg, self.body.radius_au())
    }
}

impl<'a, P> OccultGlobalEphemeris for SealedOccultGlobalServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    /// Geocentric apparent ecliptic `(lon, lat)` of the occulted body,
    /// degrees (planets via `calc_ut`, stars via `fixstar_ut`).
    fn body_lonlat_deg(&self, jd_ut: f64) -> Result<(f64, f64)> {
        self.body_lonlat_impl(jd_ut)
    }

    /// Geocentric apparent ecliptic `(lon, lat)` of the Moon, degrees.
    fn moon_lonlat_deg(&self, jd_ut: f64) -> Result<(f64, f64)> {
        self.moon_lonlat_impl(jd_ut)
    }

    /// Geocentric separation minus the sum of the apparent radii
    /// degrees.
    fn sep_minus_radii_deg(&self, jd_ut: f64) -> Result<f64> {
        self.sep_minus_radii_impl(jd_ut)
    }

    /// Shadow-axis ground geometry at `jd_ut`.
    fn where_core(&self, jd_ut: f64) -> Result<WhereGeometry> {
        self.where_geometry_at_ut(jd_ut)
    }

    /// Wrapped body-minus-Moon right-ascension gap, degrees.
    fn ra_gap_deg(&self, jd_ut: f64) -> Result<f64> {
        self.ra_gap_impl(jd_ut)
    }
}

/// Next global lunar occultation from `jd_start_ut` over the sealed
/// serving: the release `lun_occult_when_glob` `(retflag, tret)`
/// composition with the 10-slot global phase layout (`[0]` maximum,
/// `[1]` noon-transit, `[2]`/`[3]` begin/end, `[4]`/`[5]` totality
/// begin/end, `[6]`/`[7]` center-line begin/end). The body id, solar
/// character and refusal label derive from the bound leg; `ecltype`
/// is the raw public filter, `backward` selects the search direction
/// and `one_try` the single-conjunction probe (a miss returns
/// `(0, tret)` with the continuation date in `tret[0]`, exactly as in
/// the release). A walk that steps off the ephemeris is
/// [`Error::EclipseSearchExhausted`] (never a fabricated event); any
/// other provider failure propagates.
pub fn lun_occult_when_glob_sealed<P: crate::positions::StateProvider>(
    serving: &SealedOccultGlobalServing<'_, P>,
    jd_start_ut: f64,
    ecltype: i32,
    backward: bool,
    one_try: bool,
) -> Result<(i32, [f64; 10])> {
    let leg = serving.body_leg();
    lun_occult_when_glob_via(
        serving,
        jd_start_ut,
        leg.body_id_opt(),
        ecltype,
        leg.is_sun(),
        &leg.label(),
        backward,
        one_try,
    )
}

impl<'a, P> SealedOccultLocalServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    /// Checked constructor for the planet leg: borrows the sealed
    /// serving pair, the Delta-T provider and the sealed nutation
    /// source (the same source bound into `serving`), and validates
    /// the observer triple. Flag selections outside the verified
    /// tropical basis are refused (same gate as the eclipse
    /// servings). A negative body id folds to the Sun. No I/O, no
    /// globals.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_planet(
        serving: &'a crate::positions::SealedServing<'a, P>,
        star_serving: &'a crate::stars::SealedStarServing<'a, P>,
        deltat_at: &'a dyn Fn(f64) -> Result<f64>,
        nutation_at: &'a dyn Fn(f64) -> Result<(f64, f64)>,
        calc_flags: i32,
        body_id: i32,
        lon_deg_east: f64,
        lat_deg: f64,
        alt_m: f64,
    ) -> Result<Self> {
        let body_id = if body_id < 0 {
            crate::constants::SUN
        } else {
            body_id
        };
        Self::checked(
            serving,
            star_serving,
            deltat_at,
            nutation_at,
            calc_flags,
            OccultBodyLeg::Planet(body_id),
            lon_deg_east,
            lat_deg,
            alt_m,
            "occ_loc_serving",
        )
    }

    /// Checked constructor for the star leg: resolves the star name
    /// through the v1 (`fixstar_ut`) grammar (same resolution as
    /// [`SealedStarEventServing`]).
    /// No I/O, no globals.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_star(
        serving: &'a crate::positions::SealedServing<'a, P>,
        star_serving: &'a crate::stars::SealedStarServing<'a, P>,
        deltat_at: &'a dyn Fn(f64) -> Result<f64>,
        nutation_at: &'a dyn Fn(f64) -> Result<(f64, f64)>,
        calc_flags: i32,
        star_name: &str,
        lon_deg_east: f64,
        lat_deg: f64,
        alt_m: f64,
    ) -> Result<Self> {
        let star = match crate::stars::resolve_star_ref(star_name) {
            Ok(matched) => matched,
            Err(crate::stars::StarResolveError::Miss(failure)) => {
                return Err(failure.to_error(star_name));
            }
            Err(crate::stars::StarResolveError::Refused(err)) => return Err(err),
        };
        Self::checked(
            serving,
            star_serving,
            deltat_at,
            nutation_at,
            calc_flags,
            OccultBodyLeg::Star {
                star,
                name: star_name.to_string(),
            },
            lon_deg_east,
            lat_deg,
            alt_m,
            "occ_loc_serving",
        )
    }

    /// Shared checked construction (flag gate, site validation plus
    /// Delta-T boundary words).
    #[allow(clippy::too_many_arguments)]
    fn checked(
        serving: &'a crate::positions::SealedServing<'a, P>,
        star_serving: &'a crate::stars::SealedStarServing<'a, P>,
        deltat_at: &'a dyn Fn(f64) -> Result<f64>,
        nutation_at: &'a dyn Fn(f64) -> Result<(f64, f64)>,
        calc_flags: i32,
        body: OccultBodyLeg,
        lon_deg_east: f64,
        lat_deg: f64,
        alt_m: f64,
        entry: &'static str,
    ) -> Result<Self> {
        use crate::constants::{
            FLG_BARYCTR, FLG_HELCTR, FLG_J2000, FLG_NONUT, FLG_RADIANS, FLG_SIDEREAL, FLG_XYZ,
        };
        if (calc_flags & (FLG_SIDEREAL | FLG_J2000 | FLG_NONUT | FLG_XYZ | FLG_RADIANS)) != 0
            || (calc_flags & (FLG_HELCTR | FLG_BARYCTR)) != 0
        {
            return Err(Error::UnsupportedFlagCombination {
                entry,
                flags: calc_flags,
            });
        }
        validate_geopos(lon_deg_east, lat_deg, alt_m)?;
        let dt_lo = deltat_at(crate::time::SIDEREAL_MODERN_START_JD)?;
        if !dt_lo.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt_lo });
        }
        let dt_hi = deltat_at(crate::time::SIDEREAL_MODERN_END_JD)?;
        if !dt_hi.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt_hi });
        }
        let (offset_lo_deg, offset_hi_deg) =
            crate::time::sidereal_continuity_offsets(dt_lo, dt_hi)?;
        Ok(Self {
            serving,
            star_serving,
            deltat_at,
            nutation_at,
            calc_flags,
            offset_lo_deg,
            offset_hi_deg,
            body,
            lon_deg_east,
            lat_deg,
            alt_m,
        })
    }

    /// Bound body leg.
    pub fn body_leg(&self) -> &OccultBodyLeg {
        &self.body
    }

    /// Bound observer triple `(longitude_deg_east, latitude_deg,
    /// altitude_m)` in release `geopos` order.
    pub fn observer_lon_lat_alt(&self) -> (f64, f64, f64) {
        (self.lon_deg_east, self.lat_deg, self.alt_m)
    }

    /// TT epoch for a UT epoch through the bound Delta-T provider.
    pub fn jd_tt_for_ut(&self, jd_ut: f64) -> Result<f64> {
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let dt = (self.deltat_at)(jd_ut)?;
        if !dt.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt });
        }
        Ok(jd_ut + dt)
    }

    /// Map a sealed provider miss onto the public range contract
    /// a missing channel inside asset
    /// coverage stays the typed miss, past coverage it names the
    /// requested body with the asset window.
    fn escalate(&self, jd_tt: f64, body_id: i32, err: Error) -> Error {
        match err {
            Error::LebBodyNotFound { .. } | Error::LebOutOfRange { .. } => {
                crate::positions::escalate_sealed_provider_miss(
                    crate::context::ProviderMode::SealedLeb,
                    self.serving.provider(),
                    body_id,
                    jd_tt,
                    err,
                )
            }
            other => other,
        }
    }

    /// Apparent ecliptic words `(lon, lat, dist, dlon, dlat, ddist)` for
    /// one planet id at a TT epoch.
    fn ecliptic_words_at_tt(&self, body_id: i32, jd_tt: f64) -> Result<[f64; 6]> {
        use crate::constants::FLG_SPEED;
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let (pos, _) = crate::positions::serving_compute_tt(
            self.serving,
            jd_tt,
            body_id,
            self.calc_flags | FLG_SPEED,
        )
        .map_err(|e| self.escalate(jd_tt, body_id, e))?;
        Ok(pos)
    }

    /// Geocentric apparent ecliptic `(lon, lat)` of the occulted body,
    /// degrees.
    fn body_lonlat_impl(&self, jd_ut: f64) -> Result<(f64, f64)> {
        match &self.body {
            OccultBodyLeg::Planet(id) => {
                let jd_tt = self.jd_tt_for_ut(jd_ut)?;
                let words = self.ecliptic_words_at_tt(*id, jd_tt)?;
                Ok((words[0], words[1]))
            }
            OccultBodyLeg::Star { star, .. } => {
                use crate::constants::FLG_SPEED;
                let (place, _, _) = self.star_serving.fixstar_ut(star, jd_ut, FLG_SPEED, None)?;
                Ok((place[0], place[1]))
            }
        }
    }

    /// Geocentric apparent ecliptic `(lon, lat)` of the Moon, degrees
    fn moon_lonlat_impl(&self, jd_ut: f64) -> Result<(f64, f64)> {
        use crate::constants::MOON;
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let words = self.ecliptic_words_at_tt(MOON, jd_tt)?;
        Ok((words[0], words[1]))
    }

    /// Geocentric apparent equatorial Cartesian AU vector of the
    /// occulted body.
    fn body_equatorial_xyz_au(&self, jd_ut: f64) -> Result<(f64, f64, f64)> {
        use crate::constants::{FLG_EQUATORIAL, FLG_XYZ};
        match &self.body {
            OccultBodyLeg::Planet(id) => {
                let jd_tt = self.jd_tt_for_ut(jd_ut)?;
                let shape = self.calc_flags | FLG_EQUATORIAL;
                let (pos, _) =
                    crate::positions::serving_compute_tt(self.serving, jd_tt, *id, shape)
                        .map_err(|e| self.escalate(jd_tt, *id, e))?;
                let xyz = crate::positions::apply_output_flags(pos, FLG_XYZ);
                Ok((xyz[0], xyz[1], xyz[2]))
            }
            OccultBodyLeg::Star { star, .. } => {
                let (place, _, _) =
                    self.star_serving
                        .fixstar_ut(star, jd_ut, FLG_EQUATORIAL | FLG_XYZ, None)?;
                Ok((place[0], place[1], place[2]))
            }
        }
    }

    /// Topocentric apparent ecliptic `(lon, lat, dist)` of the
    /// occulted body.
    fn body_topo_ecliptic_at_ut(&self, jd_ut: f64) -> Result<(f64, f64, f64)> {
        match &self.body {
            OccultBodyLeg::Planet(id) => {
                use crate::constants::{FLG_SPEED, FLG_TOPOCTR};
                if !jd_ut.is_finite() {
                    return Err(Error::NonFiniteArgument { value: jd_ut });
                }
                let jd_tt = self.jd_tt_for_ut(jd_ut)?;
                let site = (self.lon_deg_east, self.lat_deg, self.alt_m);
                let c2t_at = |tt: f64, ut1: f64| crate::erfa::c2t06a(tt, ut1, 0.0, 0.0);
                let (pos, _) = crate::positions::serving_compute_tt_topocentric(
                    self.serving,
                    jd_tt,
                    jd_ut,
                    *id,
                    self.calc_flags | FLG_TOPOCTR | FLG_SPEED,
                    site,
                    &c2t_at,
                    false,
                    crate::positions::SiteEllipsoid::Wgs84,
                )
                .map_err(|e| self.escalate(jd_tt, *id, e))?;
                Ok((pos[0], pos[1], pos[2]))
            }
            OccultBodyLeg::Star { star, .. } => {
                use crate::constants::FLG_SPEED;
                let (place, _, _) = self.star_serving.fixstar_ut(star, jd_ut, FLG_SPEED, None)?;
                Ok((place[0], place[1], place[2]))
            }
        }
    }

    /// Topocentric apparent ecliptic `(lon, lat, dist)` of the Moon
    fn moon_topo_ecliptic_at_ut(&self, jd_ut: f64) -> Result<(f64, f64, f64)> {
        use crate::constants::{FLG_SPEED, FLG_TOPOCTR, MOON};
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let site = (self.lon_deg_east, self.lat_deg, self.alt_m);
        let c2t_at = |tt: f64, ut1: f64| crate::erfa::c2t06a(tt, ut1, 0.0, 0.0);
        let (pos, _) = crate::positions::serving_compute_tt_topocentric(
            self.serving,
            jd_tt,
            jd_ut,
            MOON,
            self.calc_flags | FLG_TOPOCTR | FLG_SPEED,
            site,
            &c2t_at,
            false,
            crate::positions::SiteEllipsoid::Wgs84,
        )
        .map_err(|e| self.escalate(jd_tt, MOON, e))?;
        Ok((pos[0], pos[1], pos[2]))
    }

    /// Topocentric apparent `(body_radius, moon_radius, separation)`
    /// in degrees.
    fn topo_radii_sep_impl(&self, jd_ut: f64, inner: bool) -> Result<(f64, f64, f64)> {
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let body_p = self.body_topo_ecliptic_at_ut(jd_ut)?;
        let moon_p = self.moon_topo_ecliptic_at_ut(jd_ut)?;
        let rmoon_au = if inner {
            ECL_RMOON_INNER_AU
        } else {
            ECL_RMOON_AU
        };
        let (rbody, rmoon) =
            occultation_radii_deg(self.body.radius_au(), body_p.2, rmoon_au, moon_p.2)?;
        let sep = vincenty_separation_deg(body_p.0, body_p.1, moon_p.0, moon_p.1)?;
        Ok((rbody, rmoon, sep))
    }

    /// `(azimuth_deg, true_altitude_deg, apparent_altitude_deg)` for a
    /// topocentric ecliptic place: shared `azalt`-equivalent composition
    /// (see [`altaz_of_topo_ecliptic_impl`]).
    fn altaz_of_topo_ecliptic(
        &self,
        topo_ecl: (f64, f64, f64),
        jd_ut: f64,
    ) -> Result<(f64, f64, f64)> {
        altaz_of_topo_ecliptic_impl(
            topo_ecl,
            jd_ut,
            self.lon_deg_east,
            self.lat_deg,
            self.alt_m,
            self.deltat_at,
            self.offset_lo_deg,
            self.offset_hi_deg,
        )
    }

    /// Local circumstances at `jd_ut`.
    fn local_how_impl(&self, jd_ut: f64) -> Result<[f64; 20]> {
        self.how_impl(jd_ut, false)
    }

    /// Center-point circumstances for the fixed-instant where
    /// composition: the uncapped disc-area ratio in `attr[2]` for planet
    /// targets in totality/annularity, plus the solar
    /// magnitude/saros tail for the Sun.
    pub fn local_how_attr_where(&self, jd_ut: f64) -> Result<[f64; 20]> {
        self.how_impl(jd_ut, true)
    }

    /// Shared how-core acquisition over the body/moon topo places with
    /// the body `azalt` triple; `where_convention` selects the release
    /// WHERE-path uncapped `attr[2]` ratio for planet targets.
    fn how_impl(&self, jd_ut: f64, where_convention: bool) -> Result<[f64; 20]> {
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let body_p = self.body_topo_ecliptic_at_ut(jd_ut)?;
        let moon_p = self.moon_topo_ecliptic_at_ut(jd_ut)?;
        let azalt = self.altaz_of_topo_ecliptic(body_p, jd_ut)?;
        let is_solar = self.body.is_sun();
        let saros = if is_solar {
            saros_info(jd_ut, EclipseKind::Solar)
        } else {
            (0.0, 0.0)
        };
        let how = sol_how_core_via(
            body_p,
            moon_p,
            self.body.radius_au(),
            is_solar,
            where_convention,
            azalt,
            self.alt_m,
            saros,
        )?;
        Ok(how.attr())
    }

    /// Shadow-axis ground geometry at a UT epoch.
    fn where_geometry_at_ut(&self, jd_ut: f64) -> Result<WhereGeometry> {
        use crate::time::{SiderealDeltaT, greenwich_apparent_sidereal_time_hours};
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let moon_xyz = {
            use crate::constants::{FLG_EQUATORIAL, FLG_XYZ, MOON};
            let jd_tt = self.jd_tt_for_ut(jd_ut)?;
            let shape = self.calc_flags | FLG_EQUATORIAL;
            let (pos, _) = crate::positions::serving_compute_tt(self.serving, jd_tt, MOON, shape)
                .map_err(|e| self.escalate(jd_tt, MOON, e))?;
            let xyz = crate::positions::apply_output_flags(pos, FLG_XYZ);
            (xyz[0], xyz[1], xyz[2])
        };
        let body_xyz = self.body_equatorial_xyz_au(jd_ut)?;
        let dt_days = (self.deltat_at)(jd_ut)?;
        if !dt_days.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt_days });
        }
        let jd_tt = jd_ut + dt_days;
        let (dpsi_rad, deps_rad) = (self.nutation_at)(jd_tt)?;
        let dpsi_deg = dpsi_rad.to_degrees();
        let eps_mean_deg = crate::positions::vondrak_mean_obliquity_deg(jd_tt)?;
        let eps_true_deg = eps_mean_deg + deps_rad.to_degrees();
        let sdt = SiderealDeltaT::try_new(dt_days, self.offset_lo_deg, self.offset_hi_deg)?;
        let gast_deg =
            greenwich_apparent_sidereal_time_hours(jd_ut, eps_true_deg, dpsi_deg, &sdt)? * 15.0;
        eclipse_where_core_via(moon_xyz, body_xyz, gast_deg, self.body.radius_au())
    }

    /// Body rise/set transits from `jd_from`.
    /// Stars travel the sealed star rise/set path (always disc
    /// center, as the release documents).
    fn body_rise_set_impl(&self, jd_from: f64) -> Result<(Option<f64>, Option<f64>)> {
        use crate::constants::{BIT_DISC_BOTTOM, CALC_RISE, CALC_SET};
        if !jd_from.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_from });
        }
        let one = |rsmi: i32| -> Result<Option<f64>> {
            let (rc, tret) = match &self.body {
                OccultBodyLeg::Planet(id) => {
                    let transit = SealedTransitServing::try_new(
                        self.serving,
                        self.deltat_at,
                        *id,
                        self.calc_flags,
                    )?;
                    rise_set_ut_serving(
                        &transit,
                        &RiseSetRequest {
                            jd_start_ut: jd_from,
                            lon_deg_east: self.lon_deg_east,
                            lat_deg: self.lat_deg,
                            height_m: self.alt_m,
                            atpress_mbar: 0.0,
                            attemp_c: 0.0,
                            horhgt_deg: 0.0,
                            rsmi,
                        },
                    )?
                }
                OccultBodyLeg::Star { name, .. } => {
                    let star_event = SealedStarEventServing::try_new(
                        self.star_serving,
                        name,
                        self.deltat_at,
                        self.nutation_at,
                    )?;
                    star_rise_set_ut_serving(
                        &star_event,
                        &StarRiseSetRequest {
                            jd_start_ut: jd_from,
                            lon_deg_east: self.lon_deg_east,
                            lat_deg: self.lat_deg,
                            height_m: self.alt_m,
                            atpress_mbar: 0.0,
                            attemp_c: 0.0,
                            horhgt_deg: 0.0,
                            rsmi,
                        },
                    )?
                }
            };
            Ok(if rc >= 0 { Some(tret[0]) } else { None })
        };
        Ok((
            one(CALC_RISE | BIT_DISC_BOTTOM)?,
            one(CALC_SET | BIT_DISC_BOTTOM)?,
        ))
    }

    /// Solar rise/set transits at `jd`.
    fn sun_rise_set_impl(&self, jd: f64) -> Result<(Option<f64>, Option<f64>)> {
        use crate::constants::{CALC_RISE, CALC_SET, SUN};
        if !jd.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd });
        }
        let sun =
            SealedTransitServing::try_new(self.serving, self.deltat_at, SUN, self.calc_flags)?;
        let one = |rsmi: i32| -> Result<Option<f64>> {
            let (rc, tret) = rise_set_ut_serving(
                &sun,
                &RiseSetRequest {
                    jd_start_ut: jd,
                    lon_deg_east: self.lon_deg_east,
                    lat_deg: self.lat_deg,
                    height_m: self.alt_m,
                    atpress_mbar: 0.0,
                    attemp_c: 0.0,
                    horhgt_deg: 0.0,
                    rsmi,
                },
            )?;
            Ok(if rc >= 0 { Some(tret[0]) } else { None })
        };
        Ok((one(CALC_RISE)?, one(CALC_SET)?))
    }
}

impl<'a, P> OccultLocalEphemeris for SealedOccultLocalServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    /// Geocentric apparent ecliptic `(lon, lat)` of the occulted body,
    /// degrees (conjunction-seek source; stars geocentric).
    fn body_lonlat_deg(&self, jd_ut: f64) -> Result<(f64, f64)> {
        self.body_lonlat_impl(jd_ut)
    }

    /// Geocentric apparent ecliptic `(lon, lat)` of the Moon, degrees.
    fn moon_lonlat_deg(&self, jd_ut: f64) -> Result<(f64, f64)> {
        self.moon_lonlat_impl(jd_ut)
    }

    /// Topocentric apparent `(body_radius, moon_radius, separation)` at
    /// `jd_ut`, degrees.
    fn topo_radii_sep_deg(&self, jd_ut: f64, inner: bool) -> Result<(f64, f64, f64)> {
        self.topo_radii_sep_impl(jd_ut, inner)
    }

    /// Local-circumstance altitude of the occulted body at `jd_ut`,
    /// degrees.
    fn body_alt_deg(&self, jd_ut: f64) -> Result<f64> {
        Ok(self.local_how_impl(jd_ut)?[6])
    }

    /// Local circumstances at the maximum.
    fn local_how_attr(&self, jd_ut: f64) -> Result<[f64; 20]> {
        self.local_how_impl(jd_ut)
    }

    /// Shadow-axis ground geometry at `jd_ut`.
    fn where_core(&self, jd_ut: f64) -> Result<WhereGeometry> {
        self.where_geometry_at_ut(jd_ut)
    }

    /// Body rise/set transits from `jd_from`: `Some(t)` carries
    /// `tr[0]` when `rc >= 0`, `None` when the body never rises/sets.
    fn body_rise_set(&self, jd_from: f64) -> Result<(Option<f64>, Option<f64>)> {
        self.body_rise_set_impl(jd_from)
    }

    /// Solar rise/set transits at `jd`: `Some(t)` carries `tr[0]` when `rc
    /// >= 0`, `None` otherwise.
    fn sun_rise_set(&self, jd: f64) -> Result<(Option<f64>, Option<f64>)> {
        self.sun_rise_set_impl(jd)
    }
}

/// Next local lunar occultation observable from the bound place, from
/// `jd_start_ut` over the sealed serving: the release
/// `(retflag, tret, attr)` local-occultation result
/// (`tret`: `[0]` maximum, `[1]` first contact, `[2]` second, `[3]`
/// third, `[4]` fourth, `[5]` body rise and `[6]` body set during the
/// occultation else `0.0`; `attr`: local circumstances at the maximum
/// with `attr[3]` = core-shadow width at the global shadow-center
/// point; `retflag`: phase bits plus `ECL_VISIBLE`, per-contact
/// visibility bits, daylight bits and the global `ECL_NONCENTRAL`
/// character). `backward` selects the search direction and `one_try`
/// the single-conjunction probe (a miss returns `(0, tret,
/// zero_attr)` with the continuation date in `tret[0]`, exactly as in
/// the release). A walk that steps off the ephemeris is
/// [`Error::EclipseSearchExhausted`] (never a fabricated event); any
/// other provider failure propagates.
pub fn lun_occult_when_loc_sealed<P: crate::positions::StateProvider>(
    serving: &SealedOccultLocalServing<'_, P>,
    jd_start_ut: f64,
    backward: bool,
    one_try: bool,
) -> Result<(i32, [f64; 10], [f64; 20])> {
    let body_id = serving.body_leg().body_id_opt();
    // The `None` body id selects the point-source contact coincidence
    // inside the walk, exactly as the release star-name arm does.
    lun_occult_when_loc_via(serving, jd_start_ut, body_id, backward, one_try)
}

/// Disc-overlap rule:
/// an occultation is in progress while the great-circle separation does
/// not exceed the sum of the occulting and occulted angular radii. This
/// is finite-source depth handling — a point coincidence test would miss
/// grazing and partial events.
pub fn planet_occult_threshold_met(
    separation_deg: f64,
    occulting_radius_deg: f64,
    target_radius_deg: f64,
) -> bool {
    separation_deg <= occulting_radius_deg + target_radius_deg
}

/// Occultation character from the minimum separation: total while the smaller disc is fully inside the larger
/// (`min_sep < |occ - target|`), partial otherwise.
pub fn classify_planet_occultation(
    min_sep_deg: f64,
    occulting_radius_deg: f64,
    target_radius_deg: f64,
) -> i32 {
    if min_sep_deg < (occulting_radius_deg - target_radius_deg).abs() {
        ECL_TOTAL
    } else {
        ECL_PARTIAL
    }
}

/// Request validation for both planetary-occultation searches (the
/// release `ValueError` arms and messages, surfaced here as
/// [`Error::InvalidEclipseSearch`], the release base-`Error`, exactly as
/// the occultation star gate maps its release `ValueError`). `occulted` is
/// `None` for a star target; `starname` carries the release `starname` /
/// `star_name` spelling for the missing-target message only — name
/// resolution itself stays star-catalog/minor-body serving-layer work. The
/// planet-ID membership rule mirrors the release `_PLANET_MAP` gate (the
/// classical Sun/Moon/planet set); see the sealed-serving note
/// below for the correspondence with [`body_uses_jpl_ephemeris`].
pub fn validate_planet_occult_request(
    occulting: i32,
    occulted: Option<i32>,
    starname: &str,
) -> Result<()> {
    if occulted.is_none() && starname.is_empty() {
        return Err(Error::InvalidEclipseSearch {
            entry: "planet_occult",
            detail:
                "Either occulted_planet ID or starname must be specified for occultation search"
                    .to_string(),
        });
    }
    if occulting == SUN {
        return Err(Error::InvalidEclipseSearch {
            entry: "planet_occult",
            detail: "Sun cannot be an occulting body - use _sol_eclipse_when_glob_pythonic"
                .to_string(),
        });
    }
    if occulting == MOON {
        return Err(Error::InvalidEclipseSearch {
            entry: "planet_occult",
            detail: "Moon cannot be an occulting body - use lun_occult_when_glob".to_string(),
        });
    }
    if !body_uses_jpl_ephemeris(occulting) {
        return Err(Error::InvalidEclipseSearch {
            entry: "planet_occult",
            detail: format!("Invalid occulting planet ID: {occulting}"),
        });
    }
    if let Some(target) = occulted {
        if !body_uses_jpl_ephemeris(target) {
            return Err(Error::InvalidEclipseSearch {
                entry: "planet_occult",
                detail: format!("Invalid occulted planet ID: {target}"),
            });
        }
        if target == occulting {
            return Err(Error::InvalidEclipseSearch {
                entry: "planet_occult",
                detail: "Occulting and occulted planets cannot be the same".to_string(),
            });
        }
    }
    Ok(())
}

/// Great-circle separation at one epoch over the walk provider.
pub fn planet_occult_separation_via(eph: &impl PlanetOccultEphemeris, jd_ut: f64) -> Result<f64> {
    let (occ_ra, occ_dec, _) = eph.occulting_position(jd_ut)?;
    let (tgt_ra, tgt_dec, _) = eph.target_position(jd_ut)?;
    planet_angular_separation_deg(occ_ra, occ_dec, tgt_ra, tgt_dec)
}

/// One-epoch occultation test:
/// `(is_occultation, separation, occulting_radius, target_radius)`.
pub fn check_planet_occultation_via(
    eph: &impl PlanetOccultEphemeris,
    jd_ut: f64,
) -> Result<(bool, f64, f64, f64)> {
    let (occ_ra, occ_dec, occ_radius) = eph.occulting_position(jd_ut)?;
    let (tgt_ra, tgt_dec, target_radius) = eph.target_position(jd_ut)?;
    let separation = planet_angular_separation_deg(occ_ra, occ_dec, tgt_ra, tgt_dec)?;
    let hit = planet_occult_threshold_met(separation, occ_radius, target_radius);
    Ok((hit, separation, occ_radius, target_radius))
}

/// Time of minimum separation over `[jd_lo, jd_hi]`.
pub fn planet_occult_minimum_via(
    eph: &impl PlanetOccultEphemeris,
    jd_lo: f64,
    jd_hi: f64,
) -> Result<f64> {
    let mut sep = |jd: f64| planet_occult_separation_via(eph, jd);
    golden_section_min(
        &mut sep,
        jd_lo,
        jd_hi,
        PLANET_OCCULT_GOLDEN_TOL_DAYS,
        PLANET_OCCULT_GOLDEN_ITERATIONS,
    )
}

/// Contact times around the maximum: `(first, second, third, fourth)` with `0.0` for unreached
/// contacts. The relative motion comes from `±1 h` samples whose RA
/// differences are scaled by `cos(mean dec)` into true angular
/// displacements (unchanged at `dec = 0`) and folded across the
/// 0/360-deg wrap; motion below `0.001` deg/day falls back to `0.5`
/// deg/day; both chords use the `sqrt(max(0, thr² − sep²))` guard.
/// Literature method pointers live at [`planet_occult_when_glob_sealed`].
pub fn planet_occult_contact_times_via(
    eph: &impl PlanetOccultEphemeris,
    jd_max: f64,
    min_sep_deg: f64,
    occulting_radius_deg: f64,
    target_radius_deg: f64,
) -> Result<(f64, f64, f64, f64)> {
    let dt = PLANET_OCCULT_SPEED_SAMPLE_DAYS;
    let (occ_ra1, occ_dec1) = eph.occulting_ra_dec_deg(jd_max - dt)?;
    let (occ_ra2, occ_dec2) = eph.occulting_ra_dec_deg(jd_max + dt)?;
    let (tgt_ra1, tgt_dec1) = eph.target_ra_dec_deg(jd_max - dt)?;
    let (tgt_ra2, tgt_dec2) = eph.target_ra_dec_deg(jd_max + dt)?;

    let mean_dec = (occ_dec1 + occ_dec2 + tgt_dec1 + tgt_dec2) / 4.0;
    let cos_dec = mean_dec.to_radians().cos();
    let mut d_ra = (occ_ra2 - occ_ra1) - (tgt_ra2 - tgt_ra1);
    d_ra = super::geometry::positive_mod(d_ra + 180.0, 360.0) - 180.0;
    d_ra *= cos_dec;
    let d_dec = (occ_dec2 - occ_dec1) - (tgt_dec2 - tgt_dec1);
    let mut relative_speed = (d_ra * d_ra + d_dec * d_dec).sqrt() / (2.0 * dt);
    if relative_speed < PLANET_OCCULT_SPEED_FLOOR_DEG_DAY {
        relative_speed = PLANET_OCCULT_SPEED_FALLBACK_DEG_DAY;
    }

    let outer_threshold = occulting_radius_deg + target_radius_deg;
    let inner_threshold = (occulting_radius_deg - target_radius_deg).abs();

    let (jd_first, jd_fourth) = if min_sep_deg < outer_threshold {
        let half_chord = (outer_threshold * outer_threshold - min_sep_deg * min_sep_deg)
            .max(0.0)
            .sqrt();
        let dt_outer = half_chord / relative_speed;
        (jd_max - dt_outer, jd_max + dt_outer)
    } else {
        (0.0, 0.0)
    };

    let (jd_second, jd_third) = if min_sep_deg < inner_threshold {
        let half_chord = (inner_threshold * inner_threshold - min_sep_deg * min_sep_deg)
            .max(0.0)
            .sqrt();
        let dt_inner = half_chord / relative_speed;
        (jd_max - dt_inner, jd_max + dt_inner)
    } else {
        (0.0, 0.0)
    };

    Ok((jd_first, jd_second, jd_third, jd_fourth))
}

/// Whether a provider failure ends the walk fail-closed. Any other
/// provider error propagates.
pub(super) fn planet_occult_range_ends_search(err: &Error) -> bool {
    matches!(
        err,
        Error::EphemerisRange { .. } | Error::LebOutOfRange { .. } | Error::SpkOutOfRange { .. }
    )
}

/// Human-readable target description for the exhaustion messages
pub(super) fn planet_occult_target_desc(occulted: Option<i32>, starname: &str) -> String {
    match occulted {
        None => starname.to_string(),
        Some(id) => format!("planet {id}"),
    }
}

/// Next planetary occultation globally: daily walk with
/// the `sep < 0.2` close-approach refinement (adopt a refined hit, jump
/// past a refined minimum that already lies ahead, else fall through to
/// the adaptive step), hit refinement over `±1.0` d with re-verification,
/// chord contacts, the total/partial rule, and the 10-slot `tret`
/// (`[max, 0, begin, end, total-begin, total-end, 0, 0, 0, 0]`).
/// `direction < 0` searches backward, exactly as the release
/// `direction` parameter. Returns `(ECL_TOTAL | ECL_PARTIAL, tret)`.
/// Literature method pointers live at [`planet_occult_when_glob_sealed`].
pub fn planet_occult_when_glob_via(
    eph: &impl PlanetOccultEphemeris,
    jd_start: f64,
    occulting: i32,
    occulted: Option<i32>,
    starname: &str,
    direction: i32,
) -> Result<(i32, [f64; 10])> {
    if !jd_start.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_start });
    }
    validate_planet_occult_request(occulting, occulted, starname)?;
    let sign = if direction < 0 { -1.0 } else { 1.0 };

    let mut jd = jd_start;
    for _ in 0..PLANET_OCCULT_MAX_ITERATIONS {
        let (mut is_occ, sep, _, _) = match check_planet_occultation_via(eph, jd) {
            Ok(found) => found,
            Err(e) if planet_occult_range_ends_search(&e) => break,
            Err(e) => return Err(e),
        };

        // Mutual planetary occultations subtend arcseconds and last
        // minutes: a sampled hit can fall between daily steps, so any
        // close approach is refined to its true minimum first.
        if !is_occ && sep < PLANET_OCCULT_CLOSE_SEP_DEG {
            let jd_min = planet_occult_minimum_via(
                eph,
                jd - PLANET_OCCULT_CLOSE_HALF_WINDOW_DAYS,
                jd + PLANET_OCCULT_CLOSE_HALF_WINDOW_DAYS,
            )?;
            // FIX-C28: the refinement check sits outside the release
            // `try` (only the main sample is guarded) — every
            // refinement failure propagates, fail-closed, instead of
            // ending the walk in exhaustion.
            let (hit, _, _, _) = check_planet_occultation_via(eph, jd_min)?;
            if hit {
                is_occ = true;
                jd = jd_min;
            } else if (jd_min - jd) * sign > 0.0 {
                // Past this conjunction: jump beyond it.
                jd = jd_min + sign * PLANET_OCCULT_CLOSE_HALF_WINDOW_DAYS;
                continue;
            }
        }

        if is_occ {
            let jd_max = planet_occult_minimum_via(
                eph,
                jd - PLANET_OCCULT_HIT_HALF_WINDOW_DAYS,
                jd + PLANET_OCCULT_HIT_HALF_WINDOW_DAYS,
            )?;
            // FIX-C28: same release `try` scope — the re-verification
            // check propagates every failure.
            let (verified, min_sep, occ_r, target_r) = check_planet_occultation_via(eph, jd_max)?;
            if verified {
                let (jd_first, jd_second, jd_third, jd_fourth) =
                    planet_occult_contact_times_via(eph, jd_max, min_sep, occ_r, target_r)?;
                let ecl_type = classify_planet_occultation(min_sep, occ_r, target_r);
                // Global 10-slot layout: maximum in
                // slot 0, slot 1 reserved, begin/end in slots 2/3,
                // totality begin/end in slots 4/5. NOTE: the local walk
                // below uses its own reordered `times` layout.
                let mut tret = [0.0; 10];
                tret[0] = jd_max;
                tret[2] = jd_first;
                tret[3] = jd_fourth;
                tret[4] = jd_second;
                tret[5] = jd_third;
                return Ok((ecl_type, tret));
            }
        }

        let step = if sep < PLANET_OCCULT_STEP_CLOSE_SEP_DEG {
            PLANET_OCCULT_STEP_CLOSE_DAYS
        } else if sep < PLANET_OCCULT_STEP_NEAR_SEP_DEG {
            PLANET_OCCULT_STEP_NEAR_DAYS
        } else if sep < PLANET_OCCULT_STEP_MID_SEP_DEG {
            PLANET_OCCULT_STEP_MID_DAYS
        } else {
            PLANET_OCCULT_STEP_FAR_DAYS
        };
        jd += sign * step;

        if jd > jd_start + PLANET_OCCULT_WINDOW_DAYS {
            break;
        }
    }

    Err(Error::EclipseSearchExhausted {
        entry: "planet_occult_when_glob",
        detail: format!(
            "No planetary occultation of {} by planet {occulting} found within {} years of JD {jd_start}",
            planet_occult_target_desc(occulted, starname),
            PLANET_OCCULT_MAX_SEARCH_YEARS,
        ),
    })
}

/// Refraction from true altitude, degrees.
pub fn planet_occult_apparent_altitude_deg(true_alt_deg: f64) -> f64 {
    if true_alt_deg > -1.0 {
        true_alt_deg
            + (1.02
                / (true_alt_deg + 10.3 / (true_alt_deg + 5.11))
                    .to_radians()
                    .tan())
                / 60.0
    } else {
        true_alt_deg
    }
}

/// Local-attr composition at the maximum: `magnitude = max(0, 1 − sep / (occ + target))`, `ratio =
/// occ / target` (`999.0` at or below the point-source radius),
/// `obscuration = magnitude`, then the azimuth / true altitude / apparent altitude /
/// separation slots with the remaining slots reserved zero. Note the
/// release `if target_radius > 0` else-arm is unreachable (stars carry
/// `0.0001` deg); nonpositive radii select the original literal arm.
pub fn planet_occult_local_attr(
    separation_deg: f64,
    occulting_radius_deg: f64,
    target_radius_deg: f64,
    target_az_deg: f64,
    target_alt_deg: f64,
) -> [f64; 20] {
    let magnitude = if target_radius_deg > 0.0 {
        (1.0 - separation_deg / (occulting_radius_deg + target_radius_deg)).max(0.0)
    } else if separation_deg < occulting_radius_deg {
        1.0
    } else {
        0.0
    };
    let ratio = if target_radius_deg > OCC_POINT_SOURCE_RADIUS_DEG {
        occulting_radius_deg / target_radius_deg
    } else {
        PLANET_OCCULT_POINT_SOURCE_RATIO
    };
    let mut attr = [0.0; 20];
    attr[0] = magnitude;
    attr[1] = ratio;
    attr[2] = magnitude;
    attr[4] = target_az_deg;
    attr[5] = target_alt_deg;
    attr[6] = planet_occult_apparent_altitude_deg(target_alt_deg);
    attr[7] = separation_deg;
    attr
}

/// Next planetary occultation visible from a site: up to 100 global re-searches, each candidate
/// checked over 10 interior samples plus the maximum against the
/// `-0.5`-deg floor on both bodies, then composed into the reordered
/// local `times` (`[max, first, second, third, fourth, 0, 0, 0, 0, 0]`)
/// and 20-slot `attr` layouts. Returns `(type, times, attr)`. The site
/// signature binds only `lat_deg`/`lon_deg`; the observer altitude
/// travels provider-side (bound `alt_m` in the sealed serving).
/// Literature method pointers live at [`planet_occult_when_loc_sealed`].
pub fn planet_occult_when_loc_via(
    eph: &impl PlanetOccultLocalEphemeris,
    jd_start: f64,
    occulting: i32,
    occulted: Option<i32>,
    starname: &str,
    lat_deg: f64,
    lon_deg: f64,
) -> Result<(i32, [f64; 10], [f64; 20])> {
    let result = planet_occult_when_loc_result_via(
        eph, jd_start, occulting, occulted, starname, lat_deg, lon_deg,
    )?;
    Ok((result.flags, result.times, result.attributes))
}

/// Local planetary search returning a typed binary64 result.
/// The geometry, contact layout and errors are shared with
/// [`planet_occult_when_loc_via`].
pub fn planet_occult_when_loc_result_via(
    eph: &impl PlanetOccultLocalEphemeris,
    jd_start: f64,
    occulting: i32,
    occulted: Option<i32>,
    starname: &str,
    lat_deg: f64,
    lon_deg: f64,
) -> Result<PlanetOccultLocalResult> {
    if !jd_start.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_start });
    }
    validate_planet_occult_request(occulting, occulted, starname)?;

    let mut current_jd = jd_start;
    for _ in 0..PLANET_OCCULT_MAX_GLOBAL_SEARCHES {
        let candidate = (|| {
            let (retflags, tret) =
                planet_occult_when_glob_via(eph, current_jd, occulting, occulted, starname, 0)?;
            let jd_max = tret[0];
            let jd_begin = tret[2];
            let jd_end = tret[3];
            let mut visible = false;
            if jd_begin > 0.0 && jd_end > 0.0 {
                for i in 0..=PLANET_OCCULT_LOCAL_SAMPLES {
                    let sample = jd_begin
                        + (jd_end - jd_begin) * i as f64 / PLANET_OCCULT_LOCAL_SAMPLES as f64;
                    // Every sample is consumed, including those after a hit.
                    if planet_occult_site_visible(eph, sample)? {
                        visible = true;
                    }
                }
            }
            if planet_occult_site_visible(eph, jd_max)? {
                visible = true;
            }
            if visible {
                let _occulting = eph.occulting_alt_az_deg(jd_max)?;
                let (target_alt, target_az) = eph.target_alt_az_deg(jd_max)?;
                let (occ_radius, target_radius, separation) = eph.local_radii_sep_deg(jd_max)?;
                let attr = planet_occult_local_attr(
                    separation,
                    occ_radius,
                    target_radius,
                    target_az,
                    target_alt,
                );
                let mut times = [0.0; 10];
                times[0] = jd_max;
                times[1] = if jd_begin > 0.0 { jd_begin } else { 0.0 };
                times[2] = if tret[4] > 0.0 { tret[4] } else { 0.0 };
                times[3] = if tret[5] > 0.0 { tret[5] } else { 0.0 };
                times[4] = if jd_end > 0.0 { jd_end } else { 0.0 };
                return Ok(Some(PlanetOccultLocalResult {
                    flags: retflags,
                    times,
                    attributes: attr,
                }));
            }
            current_jd = if jd_end > 0.0 {
                jd_end + 1.0
            } else {
                jd_max + 1.0
            };
            Ok(None)
        })();
        match candidate {
            Ok(Some(found)) => return Ok(found),
            Ok(None) => {}
            Err(error) => {
                if planet_occult_local_caught(&error) {
                    break;
                }
                return Err(error);
            }
        }
    }

    Err(Error::EclipseSearchExhausted {
        entry: "planet_occult_when_loc",
        detail: format!(
            "No planetary occultation of {} by planet {occulting} visible from lat={lat_deg}, lon={lon_deg} found within {} years of JD {jd_start}",
            planet_occult_target_desc(occulted, starname),
            PLANET_OCCULT_MAX_SEARCH_YEARS,
        ),
    })
}

/// Data, configuration and exhausted-search failures end a local candidate
/// walk. Geometric-domain and arithmetic failures propagate to the caller.
fn planet_occult_local_caught(error: &Error) -> bool {
    matches!(
        error.primary_error(),
        Error::ConfigurationMissing { .. }
            | Error::NetworkSealed { .. }
            | Error::EphemerisRange { .. }
            | Error::UnknownBody { .. }
            | Error::StarNotFound { .. }
            | Error::EclipseSearchExhausted { .. }
            | Error::LebOutOfRange { .. }
            | Error::SpkOutOfRange { .. }
    )
}

/// Both-body visibility test for the local walk.
pub(super) fn planet_occult_site_visible(
    eph: &impl PlanetOccultLocalEphemeris,
    jd_ut: f64,
) -> Result<bool> {
    let (occ_alt, _) = eph.occulting_alt_az_deg(jd_ut)?;
    let (target_alt, _) = eph.target_alt_az_deg(jd_ut)?;
    Ok(occ_alt > PLANET_OCCULT_MIN_ALT_DEG && target_alt > PLANET_OCCULT_MIN_ALT_DEG)
}

impl<'a, P> SealedPlanetOccultGlobalServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    /// Checked constructor for the planet-target leg: validates the
    /// release request (classical-set membership, Sun/Moon occulting
    /// refusals, same-body refusal) before the serving gates, then
    /// borrows the sealed serving pair and the Delta-T provider,
    /// resolving the sidereal-continuity offsets from the provider's
    /// own boundary words. Flag selections outside the verified
    /// tropical basis are refused (same gate as the occult servings).
    /// No I/O, no globals.
    pub fn try_new_planet(
        serving: &'a crate::positions::SealedServing<'a, P>,
        star_serving: &'a crate::stars::SealedStarServing<'a, P>,
        deltat_at: &'a dyn Fn(f64) -> Result<f64>,
        calc_flags: i32,
        occulting: i32,
        occulted: i32,
    ) -> Result<Self> {
        validate_planet_occult_request(occulting, Some(occulted), "")?;
        Self::checked(
            serving,
            star_serving,
            deltat_at,
            calc_flags,
            occulting,
            OccultBodyLeg::Planet(occulted),
            "planet_occ_glob_serving",
        )
    }

    /// Checked constructor for the star-target leg: validates the
    /// release request, then resolves the star name through the v1
    /// (`fixstar_ut`) grammar (same resolution as the occult
    /// servings). No I/O, no globals.
    pub fn try_new_star(
        serving: &'a crate::positions::SealedServing<'a, P>,
        star_serving: &'a crate::stars::SealedStarServing<'a, P>,
        deltat_at: &'a dyn Fn(f64) -> Result<f64>,
        calc_flags: i32,
        occulting: i32,
        star_name: &str,
    ) -> Result<Self> {
        validate_planet_occult_request(occulting, None, star_name)?;
        let star = match crate::stars::resolve_star_ref(star_name) {
            Ok(matched) => matched,
            Err(crate::stars::StarResolveError::Miss(failure)) => {
                return Err(failure.to_error(star_name));
            }
            Err(crate::stars::StarResolveError::Refused(err)) => return Err(err),
        };
        Self::checked(
            serving,
            star_serving,
            deltat_at,
            calc_flags,
            occulting,
            OccultBodyLeg::Star {
                star,
                name: star_name.to_string(),
            },
            "planet_occ_glob_serving",
        )
    }

    /// Shared checked construction (flag gate plus Delta-T boundary
    /// words).
    fn checked(
        serving: &'a crate::positions::SealedServing<'a, P>,
        star_serving: &'a crate::stars::SealedStarServing<'a, P>,
        deltat_at: &'a dyn Fn(f64) -> Result<f64>,
        calc_flags: i32,
        occulting: i32,
        target: OccultBodyLeg,
        entry: &'static str,
    ) -> Result<Self> {
        use crate::constants::{
            FLG_BARYCTR, FLG_HELCTR, FLG_J2000, FLG_NONUT, FLG_RADIANS, FLG_SIDEREAL, FLG_XYZ,
        };
        if (calc_flags & (FLG_SIDEREAL | FLG_J2000 | FLG_NONUT | FLG_XYZ | FLG_RADIANS)) != 0
            || (calc_flags & (FLG_HELCTR | FLG_BARYCTR)) != 0
        {
            return Err(Error::UnsupportedFlagCombination {
                entry,
                flags: calc_flags,
            });
        }
        let dt_lo = deltat_at(crate::time::SIDEREAL_MODERN_START_JD)?;
        if !dt_lo.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt_lo });
        }
        let dt_hi = deltat_at(crate::time::SIDEREAL_MODERN_END_JD)?;
        if !dt_hi.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt_hi });
        }
        let (offset_lo_deg, offset_hi_deg) =
            crate::time::sidereal_continuity_offsets(dt_lo, dt_hi)?;
        Ok(Self {
            serving,
            star_serving,
            deltat_at,
            calc_flags,
            offset_lo_deg,
            offset_hi_deg,
            occulting,
            target,
        })
    }

    /// Bound occulting planet id.
    pub fn occulting_id(&self) -> i32 {
        self.occulting
    }

    /// Bound occulted planet id (`None` for the star leg).
    pub fn occulted_id_opt(&self) -> Option<i32> {
        self.target.body_id_opt()
    }

    /// Bound target star name (empty for the planet leg; the release
    /// `starname`/`star_name` spelling).
    pub fn target_starname(&self) -> &str {
        match &self.target {
            OccultBodyLeg::Planet(_) => "",
            OccultBodyLeg::Star { name, .. } => name,
        }
    }

    /// Bound target leg.
    pub fn target_leg(&self) -> &OccultBodyLeg {
        &self.target
    }

    /// Bound sealed serving pair.
    pub fn serving(&self) -> &crate::positions::SealedServing<'a, P> {
        self.serving
    }

    /// Bound base calculation flags (without the per-acquisition shaping
    /// bits the walk call sites add).
    pub fn calc_flags(&self) -> i32 {
        self.calc_flags
    }

    /// TT epoch for a UT epoch through the bound Delta-T provider.
    pub fn jd_tt_for_ut(&self, jd_ut: f64) -> Result<f64> {
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let dt = (self.deltat_at)(jd_ut)?;
        if !dt.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt });
        }
        Ok(jd_ut + dt)
    }

    /// Map a sealed provider miss onto the public range contract
    /// a missing channel inside asset
    /// coverage stays the typed miss, past coverage it names the
    /// requested body with the asset window.
    pub(crate) fn escalate_miss(&self, jd_tt: f64, body_id: i32, err: Error) -> Error {
        match err {
            Error::LebBodyNotFound { .. } | Error::LebOutOfRange { .. } => {
                crate::positions::escalate_sealed_provider_miss(
                    crate::context::ProviderMode::SealedLeb,
                    self.serving.provider(),
                    body_id,
                    jd_tt,
                    err,
                )
            }
            other => other,
        }
    }

    /// Geocentric apparent equatorial `(RA deg, Dec deg, dist AU)` for
    /// one planet id at a TT epoch.
    fn planet_equatorial_at_tt(&self, body_id: i32, jd_tt: f64) -> Result<(f64, f64, f64)> {
        use crate::constants::{FLG_EQUATORIAL, FLG_SPEED};
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let (pos, _) = crate::positions::serving_compute_tt(
            self.serving,
            jd_tt,
            body_id,
            self.calc_flags | FLG_EQUATORIAL | FLG_SPEED,
        )
        .map_err(|e| self.escalate_miss(jd_tt, body_id, e))?;
        Ok((pos[0], pos[1], pos[2]))
    }

    /// Geocentric apparent equatorial `(RA, Dec)` of the occulting
    /// planet at a UT epoch, degrees.
    fn occulting_ra_dec_impl(&self, jd_ut: f64) -> Result<(f64, f64)> {
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let (ra, dec, _) = self.planet_equatorial_at_tt(self.occulting, jd_tt)?;
        Ok((ra, dec))
    }

    /// Geocentric apparent equatorial `(RA, Dec)` of the target at a UT
    /// epoch, degrees (planets via `calc_ut`, stars via `fixstar_ut`
    /// with `FLG_EQUATORIAL`, the release `_get_target_position`
    /// composition).
    fn target_ra_dec_impl(&self, jd_ut: f64) -> Result<(f64, f64)> {
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        match &self.target {
            OccultBodyLeg::Planet(id) => {
                let jd_tt = self.jd_tt_for_ut(jd_ut)?;
                let (ra, dec, _) = self.planet_equatorial_at_tt(*id, jd_tt)?;
                Ok((ra, dec))
            }
            OccultBodyLeg::Star { star, .. } => {
                use crate::constants::FLG_EQUATORIAL;
                let (place, _, _) =
                    self.star_serving
                        .fixstar_ut(star, jd_ut, FLG_EQUATORIAL, None)?;
                Ok((place[0], place[1]))
            }
        }
    }

    /// Dynamic angular radius of the occulting planet at a UT epoch,
    /// degrees.
    fn occulting_radius_impl(&self, jd_ut: f64) -> Result<f64> {
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let (_, _, dist_au) = self.planet_equatorial_at_tt(self.occulting, jd_tt)?;
        Ok(calc_planet_angular_radius_deg(self.occulting, dist_au))
    }

    /// Dynamic angular radius of the target at a UT epoch, degrees
    /// (planets via `_calc_planet_angular_radius`, stars the `0.0001`
    /// point-source spelling).
    fn target_radius_impl(&self, jd_ut: f64) -> Result<f64> {
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        match &self.target {
            OccultBodyLeg::Planet(id) => {
                let jd_tt = self.jd_tt_for_ut(jd_ut)?;
                let (_, _, dist_au) = self.planet_equatorial_at_tt(*id, jd_tt)?;
                Ok(calc_planet_angular_radius_deg(*id, dist_au))
            }
            OccultBodyLeg::Star { .. } => Ok(OCC_POINT_SOURCE_RADIUS_DEG),
        }
    }
}

impl<'a, P> PlanetOccultEphemeris for SealedPlanetOccultGlobalServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    /// Geocentric RA/Dec of the occulting planet, degrees.
    fn occulting_ra_dec_deg(&self, jd_ut: f64) -> Result<(f64, f64)> {
        self.occulting_ra_dec_impl(jd_ut)
    }

    /// Geocentric RA/Dec of the occulted body, degrees.
    fn target_ra_dec_deg(&self, jd_ut: f64) -> Result<(f64, f64)> {
        self.target_ra_dec_impl(jd_ut)
    }

    /// Angular radius of the occulting planet, degrees.
    fn occulting_radius_deg(&self, jd_ut: f64) -> Result<f64> {
        self.occulting_radius_impl(jd_ut)
    }

    /// Angular radius of the occulted body, degrees (planets dynamic,
    /// stars [`OCC_POINT_SOURCE_RADIUS_DEG`]).
    fn target_radius_deg(&self, jd_ut: f64) -> Result<f64> {
        self.target_radius_impl(jd_ut)
    }
}

/// Next global planetary occultation from `jd_start_ut` over the sealed
/// serving: the release `planet_occult_when_glob` `(retflag, tret)`
/// composition with the 10-slot global phase layout (`[0]` maximum,
/// `[1]` reserved, `[2]`/`[3]` begin/end, `[4]`/`[5]` totality
/// begin/end). The occulting/target pair and the search direction
/// derive from the bound legs; `direction < 0` searches backward. A
/// walk that steps off the ephemeris is
/// [`Error::EclipseSearchExhausted`] (never a fabricated event); any
/// other provider failure propagates. The search method follows Meeus,
/// Astronomical Algorithms, ch. 9 (angular separation) with the
/// Herald & Sinnott planetary-occultation refinement cited by the
/// release entry.
pub fn planet_occult_when_glob_sealed<P: crate::positions::StateProvider>(
    serving: &SealedPlanetOccultGlobalServing<'_, P>,
    jd_start_ut: f64,
    direction: i32,
) -> Result<(i32, [f64; 10])> {
    planet_occult_when_glob_via(
        serving,
        jd_start_ut,
        serving.occulting_id(),
        serving.occulted_id_opt(),
        serving.target_starname(),
        direction,
    )
}

impl<'a, P> SealedPlanetOccultLocalServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    /// Checked constructor for the planet-target leg: validates the
    /// release request and the observer triple, then borrows the sealed
    /// serving pair and the Delta-T provider. Flag selections
    /// outside the verified tropical basis are refused (same gate as
    /// the occult servings). No I/O, no globals.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_planet(
        serving: &'a crate::positions::SealedServing<'a, P>,
        star_serving: &'a crate::stars::SealedStarServing<'a, P>,
        deltat_at: &'a dyn Fn(f64) -> Result<f64>,
        calc_flags: i32,
        occulting: i32,
        occulted: i32,
        lon_deg_east: f64,
        lat_deg: f64,
        alt_m: f64,
    ) -> Result<Self> {
        validate_planet_occult_request(occulting, Some(occulted), "")?;
        Self::checked(
            serving,
            star_serving,
            deltat_at,
            calc_flags,
            occulting,
            OccultBodyLeg::Planet(occulted),
            lon_deg_east,
            lat_deg,
            alt_m,
            "planet_occ_loc_serving",
        )
    }

    /// Checked constructor for the star-target leg: validates the
    /// release request, then resolves the star name through the v1
    /// (`fixstar_ut`) grammar (same resolution as the occult
    /// servings). No I/O, no globals.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_star(
        serving: &'a crate::positions::SealedServing<'a, P>,
        star_serving: &'a crate::stars::SealedStarServing<'a, P>,
        deltat_at: &'a dyn Fn(f64) -> Result<f64>,
        calc_flags: i32,
        occulting: i32,
        star_name: &str,
        lon_deg_east: f64,
        lat_deg: f64,
        alt_m: f64,
    ) -> Result<Self> {
        validate_planet_occult_request(occulting, None, star_name)?;
        let star = match crate::stars::resolve_star_ref(star_name) {
            Ok(matched) => matched,
            Err(crate::stars::StarResolveError::Miss(failure)) => {
                return Err(failure.to_error(star_name));
            }
            Err(crate::stars::StarResolveError::Refused(err)) => return Err(err),
        };
        Self::checked(
            serving,
            star_serving,
            deltat_at,
            calc_flags,
            occulting,
            OccultBodyLeg::Star {
                star,
                name: star_name.to_string(),
            },
            lon_deg_east,
            lat_deg,
            alt_m,
            "planet_occ_loc_serving",
        )
    }

    /// Shared checked construction (flag gate, site validation plus
    /// Delta-T boundary words via the bound global serving).
    #[allow(clippy::too_many_arguments)]
    fn checked(
        serving: &'a crate::positions::SealedServing<'a, P>,
        star_serving: &'a crate::stars::SealedStarServing<'a, P>,
        deltat_at: &'a dyn Fn(f64) -> Result<f64>,
        calc_flags: i32,
        occulting: i32,
        target: OccultBodyLeg,
        lon_deg_east: f64,
        lat_deg: f64,
        alt_m: f64,
        entry: &'static str,
    ) -> Result<Self> {
        use crate::constants::{
            FLG_BARYCTR, FLG_HELCTR, FLG_J2000, FLG_NONUT, FLG_RADIANS, FLG_SIDEREAL, FLG_XYZ,
        };
        if (calc_flags & (FLG_SIDEREAL | FLG_J2000 | FLG_NONUT | FLG_XYZ | FLG_RADIANS)) != 0
            || (calc_flags & (FLG_HELCTR | FLG_BARYCTR)) != 0
        {
            return Err(Error::UnsupportedFlagCombination {
                entry,
                flags: calc_flags,
            });
        }
        validate_geopos(lon_deg_east, lat_deg, alt_m)?;
        // The global half re-checks the serving gates (flag gate plus
        // Delta-T boundary words) over the same inputs.
        let global = SealedPlanetOccultGlobalServing::checked(
            serving,
            star_serving,
            deltat_at,
            calc_flags,
            occulting,
            target,
            entry,
        )?;
        Ok(Self {
            global,
            lon_deg_east,
            lat_deg,
            alt_m,
        })
    }

    /// Bound occulting planet id.
    pub fn occulting_id(&self) -> i32 {
        self.global.occulting_id()
    }

    /// Bound occulted planet id (`None` for the star leg).
    pub fn occulted_id_opt(&self) -> Option<i32> {
        self.global.occulted_id_opt()
    }

    /// Bound target star name (empty for the planet leg).
    pub fn target_starname(&self) -> &str {
        self.global.target_starname()
    }

    /// Bound observer triple `(longitude_deg_east, latitude_deg,
    /// altitude_m)` in release `geopos` order.
    pub fn observer_lon_lat_alt(&self) -> (f64, f64, f64) {
        (self.lon_deg_east, self.lat_deg, self.alt_m)
    }

    /// Topocentric apparent ecliptic `(lon, lat, dist)` of one planet
    /// id.
    fn planet_topo_ecliptic_at_ut(&self, body_id: i32, jd_ut: f64) -> Result<(f64, f64, f64)> {
        use crate::constants::{FLG_SPEED, FLG_TOPOCTR};
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let jd_tt = self.global.jd_tt_for_ut(jd_ut)?;
        let site = (self.lon_deg_east, self.lat_deg, self.alt_m);
        let c2t_at = |tt: f64, ut1: f64| crate::erfa::c2t06a(tt, ut1, 0.0, 0.0);
        let (pos, _) = crate::positions::serving_compute_tt_topocentric(
            self.global.serving(),
            jd_tt,
            jd_ut,
            body_id,
            self.global.calc_flags() | FLG_TOPOCTR | FLG_SPEED,
            site,
            &c2t_at,
            false,
            crate::positions::SiteEllipsoid::Wgs84,
        )
        .map_err(|e| self.global.escalate_miss(jd_tt, body_id, e))?;
        Ok((pos[0], pos[1], pos[2]))
    }

    /// Topocentric apparent ecliptic `(lon, lat, dist)` of the target
    /// (planets via the sealed topocentric path, stars geocentric via
    /// `fixstar_ut`, whose diurnal parallax is far below the
    /// occultation tolerances — the release `_planet_occult_when_loc`
    /// LEB-branch composition).
    fn target_topo_ecliptic_at_ut(&self, jd_ut: f64) -> Result<(f64, f64, f64)> {
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        match self.global.target_leg() {
            OccultBodyLeg::Planet(id) => self.planet_topo_ecliptic_at_ut(*id, jd_ut),
            OccultBodyLeg::Star { star, .. } => {
                use crate::constants::FLG_SPEED;
                let (place, _, _) = self
                    .global
                    .star_serving
                    .fixstar_ut(star, jd_ut, FLG_SPEED, None)?;
                Ok((place[0], place[1], place[2]))
            }
        }
    }

    /// `(true_altitude_deg, azimuth_deg)` for one topocentric ecliptic
    /// place through the shared `azalt`-equivalent composition (the
    /// release `_get_body_altitude` / `_get_target_altitude` true leg:
    /// `atpress = 0` carries no refraction).
    fn alt_az_of_topo(&self, topo_ecl: (f64, f64, f64), jd_ut: f64) -> Result<(f64, f64)> {
        let (azimuth_deg, true_alt_deg, _) = altaz_of_topo_ecliptic_impl(
            topo_ecl,
            jd_ut,
            self.lon_deg_east,
            self.lat_deg,
            self.alt_m,
            self.global.deltat_at,
            self.global.offset_lo_deg,
            self.global.offset_hi_deg,
        )?;
        Ok((true_alt_deg, azimuth_deg))
    }

    /// Topocentric `(occulting_radius, target_radius, separation)` in
    /// degrees.
    fn local_radii_sep_impl(&self, jd_ut: f64) -> Result<(f64, f64, f64)> {
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let occ = self.planet_topo_ecliptic_at_ut(self.global.occulting, jd_ut)?;
        let tgt = self.target_topo_ecliptic_at_ut(jd_ut)?;
        let occ_radius = calc_planet_angular_radius_deg(self.global.occulting, occ.2);
        let target_radius = match self.global.target_leg() {
            OccultBodyLeg::Planet(id) => calc_planet_angular_radius_deg(*id, tgt.2),
            OccultBodyLeg::Star { .. } => OCC_POINT_SOURCE_RADIUS_DEG,
        };
        let sep = vincenty_separation_deg(occ.0, occ.1, tgt.0, tgt.1)?;
        Ok((occ_radius, target_radius, sep))
    }
}

impl<'a, P> PlanetOccultEphemeris for SealedPlanetOccultLocalServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    /// Geocentric RA/Dec of the occulting planet, degrees (delegates
    /// to the bound global serving: one acquisition spelling).
    fn occulting_ra_dec_deg(&self, jd_ut: f64) -> Result<(f64, f64)> {
        self.global.occulting_ra_dec_deg(jd_ut)
    }

    /// Geocentric RA/Dec of the occulted body, degrees (delegates to
    /// the bound global serving).
    fn target_ra_dec_deg(&self, jd_ut: f64) -> Result<(f64, f64)> {
        self.global.target_ra_dec_deg(jd_ut)
    }

    /// Angular radius of the occulting planet, degrees (delegates to
    /// the bound global serving).
    fn occulting_radius_deg(&self, jd_ut: f64) -> Result<f64> {
        self.global.occulting_radius_deg(jd_ut)
    }

    /// Angular radius of the occulted body, degrees (delegates to the
    /// bound global serving).
    fn target_radius_deg(&self, jd_ut: f64) -> Result<f64> {
        self.global.target_radius_deg(jd_ut)
    }
}

impl<'a, P> PlanetOccultLocalEphemeris for SealedPlanetOccultLocalServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    /// Apparent true altitude and (south-based) azimuth of the
    /// occulting planet, degrees.
    fn occulting_alt_az_deg(&self, jd_ut: f64) -> Result<(f64, f64)> {
        let topo = self.planet_topo_ecliptic_at_ut(self.global.occulting, jd_ut)?;
        self.alt_az_of_topo(topo, jd_ut)
    }

    /// Apparent true altitude and (south-based) azimuth of the
    /// occulted body, degrees.
    fn target_alt_az_deg(&self, jd_ut: f64) -> Result<(f64, f64)> {
        let topo = self.target_topo_ecliptic_at_ut(jd_ut)?;
        self.alt_az_of_topo(topo, jd_ut)
    }

    /// Topocentric angular radii and separation at one epoch, degrees
    fn local_radii_sep_deg(&self, jd_ut: f64) -> Result<(f64, f64, f64)> {
        self.local_radii_sep_impl(jd_ut)
    }
}

/// Next local planetary occultation from `jd_start_ut` over the sealed
/// serving: the release `_planet_occult_when_loc_impl`
/// `(type, times, attr)` composition with the reordered local `times`
/// layout (`[max, first, second, third, fourth, 0, 0, 0, 0, 0]`) and
/// the 20-slot `attr` layout. The occulting/target pair derives from
/// the bound legs and the observer triple travels bound in the
/// serving. A walk that steps off the ephemeris is
/// [`Error::EclipseSearchExhausted`] (never a fabricated event); any
/// other provider failure propagates. The search method follows Meeus,
/// Astronomical Algorithms, ch. 9 (angular separation) with the
/// Herald & Sinnott planetary-occultation refinement cited by the
/// release entry.
pub fn planet_occult_when_loc_sealed<P: crate::positions::StateProvider>(
    serving: &SealedPlanetOccultLocalServing<'_, P>,
    jd_start_ut: f64,
) -> Result<(i32, [f64; 10], [f64; 20])> {
    let (lon_deg_east, lat_deg, _) = serving.observer_lon_lat_alt();
    planet_occult_when_loc_via(
        serving,
        jd_start_ut,
        serving.occulting_id(),
        serving.occulted_id_opt(),
        serving.target_starname(),
        lat_deg,
        lon_deg_east,
    )
}
