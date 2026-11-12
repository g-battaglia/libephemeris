//! Solar-eclipse geometry, classification and search facade.
//!
//! Private children own local circumstance composition, contact refinement and
//! family-specific provider adapters. Reexports retain the event API paths; the
//! distinct path/search acquisition routes are deliberately not consolidated.

use super::*;

mod circumstances;
mod contacts;
#[path = "solar_eclipses/fictitious_serving.rs"]
mod fictitious_serving;
mod global_serving;
mod local_serving;
mod longitude_serving;
mod minor_serving;
mod star_serving;

pub(super) use longitude_serving::check_angular_serving_flags;
pub use longitude_serving::{cross_ut_serving, solcross_ut_serving};
pub use minor_serving::is_minor_event_body;

pub use local_serving::{
    sol_eclipse_how_details_sealed, sol_eclipse_magnitude_at_loc_sealed,
    sol_eclipse_when_loc_sealed, sol_how_sealed, sol_local_max_time_sealed,
};

pub(crate) use circumstances::sol_how_core_lazy;
pub use circumstances::{sol_eclipse_how_parts, sol_eclipse_where_parts, sol_how_core_via};
pub(super) use contacts::{
    besselian_contact_caught, contact_provider_value, details_contact_block, solve_contact_opt,
};
pub use contacts::{
    bisect_sign_change, details_contact_time, find_local_contact_time, solar_contact_c1_via,
    solar_contact_c2_via, solar_contact_c3_via, solar_contact_c4_via, validate_contact_time,
};
pub use global_serving::{
    sol_contact_c1_sealed, sol_contact_c2_sealed, sol_contact_c3_sealed, sol_contact_c4_sealed,
    sol_eclipse_central_line_sealed, sol_eclipse_duration_min_sealed,
    sol_eclipse_northern_limit_sealed, sol_eclipse_path_width_sealed,
    sol_eclipse_southern_limit_sealed, sol_eclipse_when_glob_sealed, sol_eclipse_where_sealed,
};

/// Eclipse magnitude from gamma with an explicit partial-visibility limit
///
/// Near the limit the magnitude ramps linearly to zero over
/// [`NEAR_MISS_GAMMA_MARGIN`]; at or past the limit it is `0.0`;
/// otherwise `(1 - |gamma| / limit) * ratio`, clamped to `[0, 1.5]`.
pub fn calculate_eclipse_magnitude_with_limit(
    gamma: f64,
    moon_sun_ratio: f64,
    gamma_limit_partial: f64,
) -> Result<f64> {
    for v in [gamma, moon_sun_ratio, gamma_limit_partial] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    if is_near_miss_eclipse_with_limit(gamma, gamma_limit_partial) {
        let remaining = gamma_limit_partial - gamma.abs();
        if remaining <= 0.0 {
            return Ok(0.0);
        }
        let magnitude = (remaining / NEAR_MISS_GAMMA_MARGIN) * SHALLOW_ECLIPSE_MAG_THRESHOLD;
        return Ok((magnitude * moon_sun_ratio).clamp(0.0, ECLIPSE_MAGNITUDE_MAX));
    }
    if gamma.abs() >= gamma_limit_partial {
        return Ok(0.0);
    }
    Ok(((1.0 - gamma.abs() / gamma_limit_partial) * moon_sun_ratio)
        .clamp(0.0, ECLIPSE_MAGNITUDE_MAX))
}

/// Eclipse magnitude from gamma under the default partial limit
/// ([`ECLIPSE_GAMMA_LIMIT_PARTIAL`]).
pub fn calculate_eclipse_magnitude(gamma: f64, moon_sun_ratio: f64) -> Result<f64> {
    calculate_eclipse_magnitude_with_limit(gamma, moon_sun_ratio, ECLIPSE_GAMMA_LIMIT_PARTIAL)
}

/// True for very shallow eclipses.
pub fn is_shallow_eclipse(magnitude: f64) -> bool {
    magnitude < SHALLOW_ECLIPSE_MAG_THRESHOLD
}

/// True when gamma sits within [`NEAR_MISS_GAMMA_MARGIN`] of the given
/// visibility limit.
pub fn is_near_miss_eclipse_with_limit(gamma: f64, gamma_limit: f64) -> bool {
    gamma.abs() > gamma_limit - NEAR_MISS_GAMMA_MARGIN
}

/// True when gamma sits within [`NEAR_MISS_GAMMA_MARGIN`] of the default
/// partial limit.
pub fn is_near_miss_eclipse(gamma: f64) -> bool {
    is_near_miss_eclipse_with_limit(gamma, ECLIPSE_GAMMA_LIMIT_PARTIAL)
}

/// Parse `"solar"`/`"lunar"` case-insensitively.
pub fn parse_eclipse_kind(value: &str) -> Result<EclipseKind> {
    if value.eq_ignore_ascii_case("solar") {
        Ok(EclipseKind::Solar)
    } else if value.eq_ignore_ascii_case("lunar") {
        Ok(EclipseKind::Lunar)
    } else {
        Err(Error::InvalidEclipseType {
            value: value.to_string(),
        })
    }
}

/// Round the binary64 value with ties to even
/// for the saros/inex cycle counts. Rust `round` rounds half away from
/// zero; the release `round(delta / CYCLE)` needs ties-to-even, which
/// `round_ties_even` provides over the same binary value.
pub(super) fn round_half_even(v: f64) -> f64 {
    v.round_ties_even()
}

/// Saros `(series, member)` for an eclipse maximum: the reference row minimizing the residual
/// of `jd - ref_jd` against whole saros cycles wins (strict `<`, so the
/// first row wins ties); past the 2-day window the
/// [`SAROS_NO_MATCH`] sentinel pair is returned.
pub fn saros_info(jd_eclipse: f64, kind: EclipseKind) -> (f64, f64) {
    let refs: &[(f64, i32, i32)] = match kind {
        EclipseKind::Solar => SOLAR_SAROS_REFERENCES,
        EclipseKind::Lunar => LUNAR_SAROS_REFERENCES,
    };
    let mut best_series = 0;
    let mut best_member = 0;
    let mut best_residual = f64::INFINITY;
    for &(ref_jd, ref_saros, ref_member) in refs {
        let delta = jd_eclipse - ref_jd;
        let n_cycles = round_half_even(delta / SAROS_CYCLE_DAYS);
        let residual = (delta - n_cycles * SAROS_CYCLE_DAYS).abs();
        if residual < best_residual {
            best_residual = residual;
            best_series = ref_saros;
            best_member = ref_member + n_cycles as i32;
        }
    }
    if best_residual > SAROS_MATCH_WINDOW_DAYS {
        return (SAROS_NO_MATCH, SAROS_NO_MATCH);
    }
    (best_series as f64, best_member as f64)
}

/// Saros series number for an eclipse maximum.
pub fn saros_number(jd_eclipse: f64, kind: EclipseKind) -> i64 {
    saros_info(jd_eclipse, kind).0 as i64
}

/// Inex series number for an eclipse maximum: the anchor minimizing the whole-cycle residual wins; no
/// residual gate is applied (the sparse anchor table cannot discriminate
/// arbitrary dates — the result is meaningful only for a known eclipse
/// maximum, exactly as the release documents).
pub fn inex_number(jd_eclipse: f64, kind: EclipseKind) -> i64 {
    let refs: &[(f64, i32)] = match kind {
        EclipseKind::Solar => SOLAR_INEX_REFERENCES,
        EclipseKind::Lunar => LUNAR_INEX_REFERENCES,
    };
    let (ref_jd, ref_inex) = refs[0];
    let delta_days = jd_eclipse - ref_jd;
    let cycles = round_half_even(delta_days / INEX_CYCLE_DAYS);
    let mut best_series = ref_inex;
    let mut best_residual = (delta_days - cycles * INEX_CYCLE_DAYS).abs();
    for &(test_jd, test_inex) in refs {
        let delta = jd_eclipse - test_jd;
        let test_cycles = round_half_even(delta / INEX_CYCLE_DAYS);
        let residual = (delta - test_cycles * INEX_CYCLE_DAYS).abs();
        if residual < best_residual {
            best_residual = residual;
            best_series = test_inex;
        }
    }
    best_series as i64
}

/// Shadow-axis ground geometry over injected states.
///
/// `moon_xyz_au` / `sun_xyz_au` are geocentric equatorial Cartesian
/// vectors in AU; `gast_deg` is Greenwich apparent sidereal time in
/// degrees; `occulted_radius_au` is the occulted-body disc radius in AU
/// The lunar radius is the fixed `ECL_RMOON_AU`. The
/// ellipsoid is handled by stretching `z` so the Earth becomes a sphere
/// of equatorial radius, with a second pass refining the stretch with
/// the local geodetic factor at the found latitude.
pub fn eclipse_where_core_via(
    moon_xyz_au: (f64, f64, f64),
    sun_xyz_au: (f64, f64, f64),
    gast_deg: f64,
    occulted_radius_au: f64,
) -> Result<WhereGeometry> {
    for v in [
        moon_xyz_au.0,
        moon_xyz_au.1,
        moon_xyz_au.2,
        sun_xyz_au.0,
        sun_xyz_au.1,
        sun_xyz_au.2,
        gast_deg,
        occulted_radius_au,
    ] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    let shadow = shadow_axis_surface(moon_xyz_au, sun_xyz_au, occulted_radius_au);
    Ok(shadow.with_sidereal_angle(gast_deg))
}

/// Shadow-axis kernel with the sidereal clock consumed after cone/surface
/// geometry. Input vectors are geocentric equatorial AU; the callback supplies
/// GAST degrees. Algebra/layout match the checked injected-word kernel. Only
/// its eager sidereal input guard is absent; callback errors stop at that stage.
pub(crate) fn eclipse_where_core_lazy(
    moon_xyz_au: (f64, f64, f64),
    sun_xyz_au: (f64, f64, f64),
    gast: &mut dyn FnMut() -> Result<f64>,
    occulted_radius_au: f64,
) -> Result<WhereGeometry> {
    for v in [
        moon_xyz_au.0,
        moon_xyz_au.1,
        moon_xyz_au.2,
        sun_xyz_au.0,
        sun_xyz_au.1,
        sun_xyz_au.2,
        occulted_radius_au,
    ] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    let shadow = shadow_axis_surface(moon_xyz_au, sun_xyz_au, occulted_radius_au);
    let gast_deg = gast()?;
    Ok(shadow.with_sidereal_angle(gast_deg))
}

/// Pure cone/surface words carried across the sidereal acquisition boundary.
struct ShadowAxisSurface {
    rmt: (f64, f64, f64),
    dsmt: f64,
    rsun_au: f64,
    dmoon_au: f64,
    stretch: f64,
    gx: f64,
    gy: f64,
    gz: f64,
    lat_sph: f64,
    cosf1: f64,
    cosf2: f64,
    r0: f64,
    d0: f64,
    cap_d0: f64,
    no_eclipse: bool,
    retc: i32,
}

/// Two-pass ellipsoid stretch over already-admitted geocentric AU vectors.
/// No clock is acquired here: even a non-event still reaches the clock stage.
fn shadow_axis_surface(
    moon_xyz_au: (f64, f64, f64),
    sun_xyz_au: (f64, f64, f64),
    occulted_radius_au: f64,
) -> ShadowAxisSurface {
    let rsun_au = occulted_radius_au;
    let rmoon_au = ECL_RMOON_AU;
    let dmoon_au = 2.0 * rmoon_au;
    let de = ECL_REARTH_AU;
    let f = ECL_EARTH_FLATTENING;

    let rmt = moon_xyz_au;
    let rst = sun_xyz_au;
    let dsmt = dist3(rmt, rst);

    let mut stretch = 1.0 - f;
    let (mut gx, mut gy, mut gz) = (0.0, 0.0, 0.0);
    let mut lat_sph = 0.0;
    // Carry the second pass's cone words into the ground-diameter stage.
    let (mut cosf1, mut cosf2) = (1.0, 1.0);
    let (mut r0, mut d0, mut cap_d0) = (0.0, 0.0, 0.0);
    for pass in 0..2 {
        let rm = (rmt.0, rmt.1, rmt.2 / stretch);
        let rs = (rst.0, rst.1, rst.2 / stretch);
        let dm = (rm.0.powi(2) + rm.1.powi(2) + rm.2.powi(2)).sqrt();
        let (mut ex, mut ey, mut ez) = (rm.0 - rs.0, rm.1 - rs.1, rm.2 - rs.2);
        let dsm = (ex.powi(2) + ey.powi(2) + ez.powi(2)).sqrt();
        ex /= dsm;
        ey /= dsm;
        ez /= dsm;
        // Core and penumbra cone half angles use the radius difference/sum.
        let sinf1 = (rsun_au - rmoon_au) / dsm;
        cosf1 = (1.0 - sinf1 * sinf1).sqrt();
        let sinf2 = (rsun_au + rmoon_au) / dsm;
        cosf2 = (1.0 - sinf2 * sinf2).sqrt();
        let s0 = -(rm.0 * ex + rm.1 * ey + rm.2 * ez);
        let r0_sq = dm * dm - s0 * s0;
        r0 = if r0_sq > 0.0 { r0_sq.sqrt() } else { 0.0 };
        d0 = (s0 / dsm * (2.0 * rsun_au - dmoon_au) - dmoon_au) / cosf1;
        cap_d0 = (s0 / dsm * (2.0 * rsun_au + dmoon_au) + dmoon_au) / cosf2;
        // Stop at closest approach when the shadow axis misses the globe.
        let disc = s0 * s0 + de * de - dm * dm;
        let s_walk = s0 - (if disc > 0.0 { disc.sqrt() } else { 0.0 });
        gx = rm.0 + s_walk * ex;
        gy = rm.1 + s_walk * ey;
        gz = rm.2 + s_walk * ez;
        lat_sph = gz.atan2(gx.hypot(gy));
        if pass == 0 {
            // Refine the stretch for the ground point's geodetic latitude.
            let (cl, sl) = (lat_sph.cos(), lat_sph.sin());
            let one_f2 = (1.0 - f).powi(2);
            stretch = one_f2 / (cl * cl + one_f2 * sl * sl).sqrt();
        }
    }

    let mut no_eclipse = false;
    let mut retc = 0;
    if de * cosf1 >= r0 {
        retc |= ECL_CENTRAL;
    } else if r0 <= de * cosf1 + d0.abs() / 2.0 {
        retc |= ECL_NONCENTRAL;
    } else if r0 <= de * cosf2 + cap_d0 / 2.0 {
        retc |= ECL_PARTIAL | ECL_NONCENTRAL;
    } else {
        no_eclipse = true;
    }
    ShadowAxisSurface {
        rmt,
        dsmt,
        rsun_au,
        dmoon_au,
        stretch,
        gx,
        gy,
        gz,
        lat_sph,
        cosf1,
        cosf2,
        r0,
        d0,
        cap_d0,
        no_eclipse,
        retc,
    }
}

impl ShadowAxisSurface {
    /// Rotate into geographic degrees, then finish true-frame shadow widths.
    fn with_sidereal_angle(self, gast_deg: f64) -> WhereGeometry {
        let Self {
            rmt,
            dsmt,
            rsun_au,
            dmoon_au,
            stretch,
            gx,
            gy,
            gz,
            lat_sph,
            cosf1,
            cosf2,
            r0,
            d0,
            cap_d0,
            no_eclipse,
            mut retc,
        } = self;
        let mut lon_out = (gy.atan2(gx).to_degrees() - gast_deg).rem_euclid(360.0);
        if lon_out > 180.0 {
            lon_out -= 360.0;
        }
        let lat_out = lat_sph.to_degrees();

        // Shadow diameters at the ground point, from true-frame distances.
        let g_true = (gx, gy, gz * stretch);
        let s_mg = dist3(rmt, g_true);
        let dcore0 = (s_mg / dsmt * (2.0 * rsun_au - dmoon_au) - dmoon_au) * cosf1 * ECL_AU_KM;
        let dcore1 = (s_mg / dsmt * (2.0 * rsun_au + dmoon_au) + dmoon_au) * cosf2 * ECL_AU_KM;
        if !no_eclipse && (retc & ECL_PARTIAL) == 0 {
            retc |= if dcore0 > 0.0 { ECL_ANNULAR } else { ECL_TOTAL };
        }

        WhereGeometry {
            retc,
            lon_deg: lon_out,
            lat_deg: lat_out,
            dcore: [
                dcore0,
                dcore1,
                r0 * ECL_AU_KM,
                d0 * ECL_AU_KM,
                cap_d0 * ECL_AU_KM,
                cosf1,
                cosf2,
                0.0,
                0.0,
                0.0,
            ],
        }
    }
}

/// True when `err` signals the end of ephemeris coverage. Lunation walks convert this into
/// [`Error::EclipseSearchExhausted`]; every other provider failure
/// propagates (or degrades to `0.0` at the documented contact-solve
/// sites only).
pub fn is_ephemeris_boundary(err: &Error) -> bool {
    matches!(
        err.primary_error(),
        Error::EphemerisRange { .. } | Error::LebOutOfRange { .. } | Error::SpkOutOfRange { .. }
    )
}

/// Parse a search-direction name.
pub fn parse_search_direction(name: &str) -> Result<EclipseSearchDirection> {
    match name {
        "forward" => Ok(EclipseSearchDirection::Forward),
        "backward" => Ok(EclipseSearchDirection::Backward),
        "bidirectional" => Ok(EclipseSearchDirection::Bidirectional),
        _ => Err(Error::InvalidEclipseSearch {
            entry: "sol_eclipse_when_glob",
            detail: format!(
                "search_direction must be one of ('forward', 'backward', 'bidirectional'), got {name:?}"
            ),
        }),
    }
}

/// Coerce a public `backwards` argument to a plain bool.
pub fn coerce_backwards(value: BackwardsArg) -> Result<bool> {
    match value {
        BackwardsArg::Bool(b) => Ok(b),
        BackwardsArg::Int(i) => Ok(i != 0),
        BackwardsArg::Text(text) => {
            let normalized = text.trim().to_lowercase();
            match normalized.as_str() {
                "backward" | "backwards" | "true" | "1" => Ok(true),
                "forward" | "forwards" | "false" | "0" | "" => Ok(false),
                _ => Err(Error::InvalidEclipseSearch {
                    entry: "sol_eclipse_when_glob",
                    detail: format!(
                        "invalid backwards value {text:?}: expected a bool, an int, \
                         or one of 'backward'/'backwards'/'true'/'1' (backward search) \
                         or 'forward'/'forwards'/'false'/'0'/'' (forward search)"
                    ),
                }),
            }
        }
    }
}

/// Drop bit `0x8000` from eclipse-search calculation flags: the bit value is shared between
/// `ECL_ONE_TRY` and `FLG_TOPOCTR`, and in the eclipse `flags` argument
/// it can only mean the former.
pub fn strip_one_try_bit(flags: i32) -> i32 {
    flags & !ECL_ONE_TRY
}

/// Reference `ifltype` acceptance test for `sol_eclipse_when_glob`
/// mask `0` accepts any eclipse;
/// a geometry-only mask never matches; the eclipse must carry one of the
/// requested type bits; with no geometry bit only a single-type mask
/// matches; with a geometry bit the eclipse centrality must match.
pub fn sol_glob_accepts(mask: i32, retflag: i32) -> bool {
    if mask == 0 {
        return true;
    }
    let geom = mask & (ECL_CENTRAL | ECL_NONCENTRAL);
    let types = mask & (ECL_TOTAL | ECL_ANNULAR | ECL_PARTIAL | ECL_ANNULAR_TOTAL);
    if types == 0 {
        return false;
    }
    if retflag & types == 0 {
        return false;
    }
    if geom == 0 {
        return types.count_ones() == 1;
    }
    retflag & geom != 0
}

/// Reject `sol_eclipse_when_glob` type masks with no realisable event
/// central partial and non-central hybrid geometries cannot occur. Only
/// the masked valid bits participate in the comparison.
pub fn sol_glob_reject_impossible(mask: i32) -> Result<()> {
    let valid =
        ECL_CENTRAL | ECL_NONCENTRAL | ECL_TOTAL | ECL_ANNULAR | ECL_PARTIAL | ECL_ANNULAR_TOTAL;
    let m = mask & valid;
    if m == (ECL_CENTRAL | ECL_PARTIAL) {
        return Err(Error::InvalidEclipseSearch {
            entry: "sol_eclipse_when_glob",
            detail: "central partial eclipses do not exist".to_string(),
        });
    }
    if m == (ECL_NONCENTRAL | ECL_ANNULAR_TOTAL) {
        return Err(Error::InvalidEclipseSearch {
            entry: "sol_eclipse_when_glob",
            detail: "non-central hybrid (annular-total) eclipses do not exist".to_string(),
        });
    }
    Ok(())
}

impl SunMoonState {
    /// Bundle the eight consumed words; every word must be finite
    /// (providers own the finite-or-`Err` contract, as in the solver cores).
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        sun_lon_deg: f64,
        sun_lat_deg: f64,
        sun_dist_au: f64,
        sun_speed_deg_day: f64,
        moon_lon_deg: f64,
        moon_lat_deg: f64,
        moon_dist_au: f64,
        moon_speed_deg_day: f64,
    ) -> Result<Self> {
        for value in [
            sun_lon_deg,
            sun_lat_deg,
            sun_dist_au,
            sun_speed_deg_day,
            moon_lon_deg,
            moon_lat_deg,
            moon_dist_au,
            moon_speed_deg_day,
        ] {
            if !value.is_finite() {
                return Err(Error::NonFiniteArgument { value });
            }
        }
        Ok(Self {
            sun_lon_deg,
            sun_lat_deg,
            sun_dist_au,
            sun_speed_deg_day,
            moon_lon_deg,
            moon_lat_deg,
            moon_dist_au,
            moon_speed_deg_day,
        })
    }
    /// Solar apparent ecliptic longitude, degrees.
    pub fn sun_lon_deg(&self) -> f64 {
        self.sun_lon_deg
    }
    /// Solar apparent ecliptic latitude, degrees.
    pub fn sun_lat_deg(&self) -> f64 {
        self.sun_lat_deg
    }
    /// Solar apparent distance, AU.
    pub fn sun_dist_au(&self) -> f64 {
        self.sun_dist_au
    }
    /// Solar ecliptic longitude speed, deg/day.
    pub fn sun_speed_deg_day(&self) -> f64 {
        self.sun_speed_deg_day
    }
    /// Lunar apparent ecliptic longitude, degrees.
    pub fn moon_lon_deg(&self) -> f64 {
        self.moon_lon_deg
    }
    /// Lunar apparent ecliptic latitude, degrees.
    pub fn moon_lat_deg(&self) -> f64 {
        self.moon_lat_deg
    }
    /// Lunar apparent distance, AU.
    pub fn moon_dist_au(&self) -> f64 {
        self.moon_dist_au
    }
    /// Lunar ecliptic longitude speed, deg/day.
    pub fn moon_speed_deg_day(&self) -> f64 {
        self.moon_speed_deg_day
    }
}

/// Fold an angle to `(-180, 180]` degrees.
pub(super) fn fold_180_deg(angle_deg: f64) -> f64 {
    let mut wrapped = angle_deg.rem_euclid(360.0);
    if wrapped > 180.0 {
        wrapped -= 360.0;
    }
    wrapped
}

/// Next New Moon (Sun-Moon conjunction) after `jd_start`.
pub fn find_next_new_moon(
    state: &mut impl FnMut(f64) -> Result<SunMoonState>,
    jd_start: f64,
) -> Result<f64> {
    find_next_syzygy(state, jd_start, 0.0)
}

/// Previous New Moon before `jd_start`.
pub fn find_previous_new_moon(
    state: &mut impl FnMut(f64) -> Result<SunMoonState>,
    jd_start: f64,
) -> Result<f64> {
    find_previous_syzygy(state, jd_start, 0.0)
}

/// Next Full Moon (Sun-Moon opposition) after `jd_start`.
pub fn find_next_full_moon(
    state: &mut impl FnMut(f64) -> Result<SunMoonState>,
    jd_start: f64,
) -> Result<f64> {
    find_next_syzygy(state, jd_start, 180.0)
}

/// Previous Full Moon before `jd_start`.
pub fn find_previous_full_moon(
    state: &mut impl FnMut(f64) -> Result<SunMoonState>,
    jd_start: f64,
) -> Result<f64> {
    find_previous_syzygy(state, jd_start, 180.0)
}

/// Angular distance of the Moon from the nearest lunar node, degrees
pub fn moon_node_distance_deg(node_lon_deg: f64, moon_lon_deg: f64) -> Result<f64> {
    if !node_lon_deg.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: node_lon_deg,
        });
    }
    if !moon_lon_deg.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: moon_lon_deg,
        });
    }
    let dist_north = (fold_180_deg(moon_lon_deg - node_lon_deg)).abs();
    let dist_south = (fold_180_deg(moon_lon_deg - (node_lon_deg + 180.0))).abs();
    Ok(dist_north.min(dist_south))
}

/// Golden-section minimum of `objective` on `[lo, hi]`. Provider errors
/// propagate: a minimum over a coverage gap must never fabricate.
pub fn golden_section_min(
    objective: &mut impl FnMut(f64) -> Result<f64>,
    lo: f64,
    hi: f64,
    tolerance_days: f64,
    max_iter: usize,
) -> Result<f64> {
    // Single-loop spelling: the traced core below runs this exact
    // sequence with a discarded sink, so every caller observes the
    // historical words unchanged.
    let mut sink = Vec::new();
    golden_section_min_traced(objective, lo, hi, tolerance_days, max_iter, &mut sink)
}

/// Golden-section minimization with a per-decision word trace
/// (same contract as [`golden_section_min`): the loop body is shared —
/// [`golden_section_min`] delegates here with a discarded sink — so the
/// returned minimum is bitwise identical to the untraced call on the
/// same objective words. Each pushed [`GoldenTraceRow`] records the
/// decision inputs before the arm update; the trace length is the
/// number of executed decisions (at most `max_iter`). The sink push
/// cannot alter the floating-point words: it runs after the arm
/// comparison on already-computed values.
pub fn golden_section_min_traced(
    objective: &mut impl FnMut(f64) -> Result<f64>,
    lo: f64,
    hi: f64,
    tolerance_days: f64,
    max_iter: usize,
    trace: &mut Vec<GoldenTraceRow>,
) -> Result<f64> {
    if !lo.is_finite() {
        return Err(Error::NonFiniteArgument { value: lo });
    }
    if !hi.is_finite() {
        return Err(Error::NonFiniteArgument { value: hi });
    }
    if hi <= lo {
        return Err(Error::InvalidSearchInterval { a: lo, b: hi });
    }
    golden_section_min_source(objective, lo, hi, tolerance_days, max_iter, trace)
}

/// The consumed binary64 golden loop without bracket admission. Reversed,
/// zero-width and IEEE nonfinite bounds follow comparisons and provider reads;
/// the checked public minimizer retains its guards before calling this loop.
/// Objective errors propagate at the actual sample, before any later sample.
pub(crate) fn golden_section_min_source(
    objective: &mut impl FnMut(f64) -> Result<f64>,
    lo: f64,
    hi: f64,
    tolerance_days: f64,
    max_iter: usize,
    trace: &mut Vec<GoldenTraceRow>,
) -> Result<f64> {
    let phi = (1.0 + 5.0f64.sqrt()) / 2.0;
    let (mut jd_low, mut jd_high) = (lo, hi);
    let mut jd_a = jd_high - (jd_high - jd_low) / phi;
    let mut jd_b = jd_low + (jd_high - jd_low) / phi;
    let mut f_a = objective(jd_a)?;
    let mut f_b = objective(jd_b)?;
    for level in 0..max_iter {
        let arm_a_taken = f_a < f_b;
        let (row_a, row_b, row_fa, row_fb) = (jd_a, jd_b, f_a, f_b);
        if arm_a_taken {
            jd_high = jd_b;
            jd_b = jd_a;
            f_b = f_a;
            jd_a = jd_high - (jd_high - jd_low) / phi;
            f_a = objective(jd_a)?;
        } else {
            jd_low = jd_a;
            jd_a = jd_b;
            f_a = f_b;
            jd_b = jd_low + (jd_high - jd_low) / phi;
            f_b = objective(jd_b)?;
        }
        trace.push(GoldenTraceRow {
            level,
            jd_a: row_a,
            jd_b: row_b,
            f_a: row_fa,
            f_b: row_fb,
            arm_a_taken,
            width_after: jd_high - jd_low,
        });
        if jd_high - jd_low < tolerance_days {
            break;
        }
    }
    Ok((jd_low + jd_high) / 2.0)
}

/// Three-term Neumaier sum: retain the addition residual in a separate
/// low word, then round `high + low`. Used for fundamental-plane dot
/// products whose cancellation is sensitive to ordinary left-to-right
/// summation. Products themselves are rounded binary64 operations;
/// this is not an exact dot product or an astronomical error bound.
pub(super) fn compensated_sum3(terms: [f64; 3]) -> f64 {
    let mut high: f64 = 0.0;
    let mut low = 0.0;
    for x in terms {
        let next = high + x;
        low += if high.abs() >= x.abs() {
            (high - next) + x
        } else {
            (x - next) + high
        };
        high = next;
    }
    if low != 0.0 && low.is_finite() {
        high + low
    } else {
        high
    }
}

/// Penumbral shadow radius `l1` on the fundamental plane, Earth radii
pub fn penumbra_limit_er(sun_dist_au: f64, moon_dist_au: f64) -> Result<f64> {
    if !sun_dist_au.is_finite() {
        return Err(Error::NonFiniteArgument { value: sun_dist_au });
    }
    if !moon_dist_au.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: moon_dist_au,
        });
    }
    let au_km = ECL_AU_KM;
    let sun_dist_km = sun_dist_au * au_km;
    let moon_dist_km = moon_dist_au * au_km;
    let sun_angular_rad = ECLIPSE_SUN_RADIUS_KM / sun_dist_km;
    let moon_angular_rad = ECLIPSE_MOON_RADIUS_KM / moon_dist_km;
    let f1_rad = sun_angular_rad + moon_angular_rad;
    let l1_km = moon_dist_km * f1_rad;
    Ok(l1_km / ECLIPSE_EARTH_RADIUS_KM)
}

/// Umbral/antumbral shadow radius `l2` on the fundamental plane, Earth
/// radii: negative for umbra (total), positive for antumbra
/// (annular). The cone angle uses the Sun-Moon separation, matching the
/// sibling Besselian models.
pub fn umbra_limit_er(sun_dist_au: f64, moon_dist_au: f64) -> Result<f64> {
    if !sun_dist_au.is_finite() {
        return Err(Error::NonFiniteArgument { value: sun_dist_au });
    }
    if !moon_dist_au.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: moon_dist_au,
        });
    }
    let au_km = ECL_AU_KM;
    let sun_dist_km = sun_dist_au * au_km;
    let moon_dist_km = moon_dist_au * au_km;
    let sun_moon_dist_km = sun_dist_km - moon_dist_km;
    let f2 = ((ECLIPSE_SUN_RADIUS_KM - ECLIPSE_MOON_RADIUS_KM) / sun_moon_dist_km).atan();
    let umbra_vertex_km = if f2 > 1e-10 {
        ECLIPSE_MOON_RADIUS_KM / f2.tan()
    } else {
        f64::INFINITY
    };
    if moon_dist_km < umbra_vertex_km {
        let l2_km = ECLIPSE_MOON_RADIUS_KM - moon_dist_km * f2.tan();
        Ok(-l2_km / ECLIPSE_EARTH_RADIUS_KM)
    } else {
        let l2_km = moon_dist_km * f2.tan() - ECLIPSE_MOON_RADIUS_KM;
        Ok(l2_km / ECLIPSE_EARTH_RADIUS_KM)
    }
}

/// Refined solar-eclipse maximum: golden-section minimum of gamma over
/// `jd_approx ± range`. Returns the minimum instant; the caller re-reads
/// gamma there when it needs the value (as the release does).
pub fn refine_solar_eclipse_maximum(
    gamma: &mut impl FnMut(f64) -> Result<f64>,
    jd_approx: f64,
    search_range_days: f64,
) -> Result<f64> {
    golden_section_min(
        gamma,
        jd_approx - search_range_days,
        jd_approx + search_range_days,
        SOLAR_MAXIMUM_TOLERANCE_DAYS,
        GOLDEN_SECTION_ITERATIONS,
    )
}

/// Physical hybrid (annular-total) test: a
/// hybrid eclipse is total around maximum and annular near the path
/// ends, so the ground-point core-shadow width changes sign between
/// the umbral contacts. Outside the `[-0.001, 0.006]` Earth-radii band
/// around zero `l2` no sign change is possible (fast reject). The
/// umbral contacts come from the legacy gamma spelling
/// ([`find_contact_time_gamma`]); the widths are the ground-point core
/// diameters (`dcore[0]` words). `ground_width` read failures skip the
/// sample, except coverage misses, which propagate.
pub fn is_hybrid_solar_eclipse(
    jd_max: f64,
    l2_er: f64,
    gamma_contact: &mut impl FnMut(f64, f64, bool) -> Result<f64>,
    ground_width_km: &mut impl FnMut(f64) -> Result<f64>,
) -> Result<bool> {
    if !jd_max.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_max });
    }
    if !l2_er.is_finite() {
        return Err(Error::NonFiniteArgument { value: l2_er });
    }
    if l2_er < HYBRID_L2_FAST_REJECT_NEG || l2_er > HYBRID_L2_FAST_REJECT_POS {
        return Ok(false);
    }
    let umbral_limit = 1.0 + l2_er.abs();
    let t_u1 = gamma_contact(jd_max, umbral_limit, true)?;
    let t_u4 = gamma_contact(jd_max, umbral_limit, false)?;
    let sample_times: Vec<f64> = if t_u1 != 0.0 && t_u4 != 0.0 && t_u4 > t_u1 {
        [0.0, 0.25, 0.5, 0.75, 1.0]
            .into_iter()
            .map(|f| t_u1 + f * (t_u4 - t_u1))
            .collect()
    } else {
        vec![jd_max]
    };
    let mut widths: Vec<f64> = Vec::new();
    for t in sample_times {
        match ground_width_km(t) {
            Ok(width) => widths.push(width),
            Err(err) => {
                if is_ephemeris_boundary(&err) {
                    return Err(err);
                }
            }
        }
    }
    Ok(!widths.is_empty()
        && widths.iter().cloned().fold(f64::INFINITY, f64::min) < 0.0
        && widths.iter().cloned().fold(f64::NEG_INFINITY, f64::max) > 0.0)
}

/// Classify a solar eclipse at its refined maximum from the consumed
/// ephemeris words. `state` carries the apparent Sun/Moon
/// `(lon, lat, dist)`; `where_gate` provides the where-core retc and
/// the where-magnitude (consulted only on the near-miss arm, as in the
/// release); `hybrid_probe` answers the physical hybrid test
/// (consulted only when the core shadow reaches a central Earth).
/// Provider errors propagate.
#[allow(clippy::too_many_arguments)]
pub fn classify_solar_eclipse(
    state: &SunMoonState,
    gamma: f64,
    l1_er: f64,
    l2_er: f64,
    where_gate: &mut impl FnMut() -> Result<(i32, f64)>,
    hybrid_probe: &mut impl FnMut() -> Result<bool>,
) -> Result<SolarEclipseClassification> {
    for value in [gamma, l1_er, l2_er] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    // This provider supplies both where roles together. Adapt only that
    // boundary; the classifier still acquires it at the selected branch.
    let mut where_role = |_| where_gate();
    classify_solar_eclipse_stages(state, gamma, l1_er, l2_er, &mut where_role, hybrid_probe)
}

/// Classify the eclipse (if any) at one New Moon: the coarse node
/// pre-filter, Besselian maximum refinement, classification at the
/// refined maximum and the `mask` acceptance test. Shared by the
/// pre-probes and the main walk.
pub(super) fn check_new_moon_for_eclipse<E: SolarGlobalEphemeris + ?Sized>(
    eph: &E,
    jd_new_moon: f64,
    mask: i32,
) -> Result<Option<(i32, [f64; 10])>> {
    let moon_lon = eph.candidate_moon_longitude_deg(jd_new_moon)?;
    let node_lon = eph.mean_node_lon_deg(jd_new_moon)?;
    let node_dist = moon_node_distance_deg(node_lon, moon_lon)?;
    if node_dist < SOLAR_NODE_PREFILTER_DEG {
        let jd_max = {
            let mut gamma_at = |jd: f64| eph.gamma_er(jd);
            refine_solar_eclipse_maximum(
                &mut gamma_at,
                jd_new_moon,
                SOLAR_MAXIMUM_SEARCH_RANGE_DAYS,
            )?
        };
        let state_max = eph.sun_moon_state(jd_max)?;
        let gamma_max = eph.gamma_er(jd_max)?;
        let (l1, l2) = eph.shadow_limits_er(jd_max)?;
        let mut where_at_max = || eph.where_gate(jd_max);
        let mut hybrid_at_max = || {
            let mut contact = |jd: f64, target: f64, before: bool| -> Result<f64> {
                let mut gamma_at = |jd: f64| eph.gamma_er(jd);
                find_contact_time_gamma(
                    &mut gamma_at,
                    jd,
                    target,
                    before,
                    HYBRID_CONTACT_SEARCH_RANGE_DAYS,
                )
            };
            let mut width = |jd: f64| eph.ground_core_width_km(jd);
            is_hybrid_solar_eclipse(jd_max, l2, &mut contact, &mut width)
        };
        let classified = classify_solar_eclipse(
            &state_max,
            gamma_max,
            l1,
            l2,
            &mut where_at_max,
            &mut hybrid_at_max,
        )?;
        if classified.eclipse_type != 0 && sol_glob_accepts(mask, classified.eclipse_type) {
            let mut residuals = |jd: f64| eph.contact_residuals_km(jd);
            let mut x = |jd: f64| eph.besselian_x_er(jd);
            let times = calculate_solar_phases_besselian(
                jd_max,
                classified.eclipse_type,
                &mut residuals,
                &mut x,
            )?;
            return Ok(Some((classified.eclipse_type, times)));
        }
    }
    Ok(None)
}

/// Next global solar eclipse from `jd_start`: lunation-by-lunation walk with the
/// coarse node pre-filter, Besselian maximum refinement,
/// classification at the refined maximum, the `mask` acceptance test
/// (`0` accepts any type), the forward/backward/bidirectional
/// pre-probes and the epoch-margin direction invariants. Returns
/// `(retflag, tret)` with the 10-slot global phase layout. A walk that
/// steps off the ephemeris is [`Error::EclipseSearchExhausted`]
/// (never a fabricated event); any other provider failure propagates.
pub fn sol_eclipse_when_glob_via(
    eph: &impl SolarGlobalEphemeris,
    jd_start: f64,
    mask: i32,
    direction: EclipseSearchDirection,
) -> Result<(i32, [f64; 10])> {
    if !jd_start.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_start });
    }
    sol_glob_reject_impossible(mask)?;

    if direction == EclipseSearchDirection::Bidirectional {
        let mut state = |jd: f64| eph.sun_moon_state(jd);
        match find_previous_new_moon(&mut state, jd_start) {
            Ok(jd_prev) => {
                if jd_start - jd_prev <= ECLIPSE_PRE_PROBE_WINDOW_DAYS {
                    match check_new_moon_for_eclipse(eph, jd_prev, mask) {
                        Ok(Some(result)) => {
                            if result.1[0] > jd_start + ECLIPSE_WHEN_EPOCH_MARGIN_DAYS {
                                return Ok(result);
                            }
                        }
                        Ok(None) => {}
                        Err(err) => {
                            if !is_ephemeris_boundary(&err) {
                                return Err(err);
                            }
                        }
                    }
                }
            }
            Err(err) => {
                if !is_ephemeris_boundary(&err) {
                    return Err(err);
                }
            }
        }
    }

    if direction == EclipseSearchDirection::Backward {
        let mut state = |jd: f64| eph.sun_moon_state(jd);
        match find_next_new_moon(&mut state, jd_start) {
            Ok(jd_next) => {
                if jd_next - jd_start <= ECLIPSE_PRE_PROBE_WINDOW_DAYS {
                    match check_new_moon_for_eclipse(eph, jd_next, mask) {
                        Ok(Some(result)) => {
                            if result.1[0] < jd_start - ECLIPSE_WHEN_EPOCH_MARGIN_DAYS {
                                return Ok(result);
                            }
                        }
                        Ok(None) => {}
                        Err(err) => {
                            if !is_ephemeris_boundary(&err) {
                                return Err(err);
                            }
                        }
                    }
                }
            }
            Err(err) => {
                if !is_ephemeris_boundary(&err) {
                    return Err(err);
                }
            }
        }
        let mut jd = jd_start;
        for _ in 0..ECLIPSE_WHEN_MAX_LUNATIONS {
            let jd_prev = {
                let mut state = |jd: f64| eph.sun_moon_state(jd);
                match find_previous_new_moon(&mut state, jd) {
                    Ok(jd_prev) => jd_prev,
                    Err(err) => {
                        if is_ephemeris_boundary(&err) {
                            break;
                        }
                        return Err(err);
                    }
                }
            };
            match check_new_moon_for_eclipse(eph, jd_prev, mask) {
                Ok(Some(result)) => {
                    if result.1[0] < jd_start - ECLIPSE_WHEN_EPOCH_MARGIN_DAYS {
                        return Ok(result);
                    }
                }
                Ok(None) => {}
                Err(err) => {
                    if is_ephemeris_boundary(&err) {
                        break;
                    }
                    return Err(err);
                }
            }
            jd = jd_prev - 1.0;
        }
        return Err(Error::EclipseSearchExhausted {
            entry: "sol_eclipse_when_glob",
            detail: format!(
                "No matching solar eclipse found before the ephemeris boundary searching back from JD {jd_start}"
            ),
        });
    }

    if direction == EclipseSearchDirection::Forward {
        let mut state = |jd: f64| eph.sun_moon_state(jd);
        match find_previous_new_moon(&mut state, jd_start) {
            Ok(jd_prev) => {
                if jd_start - jd_prev <= ECLIPSE_PRE_PROBE_WINDOW_DAYS {
                    match check_new_moon_for_eclipse(eph, jd_prev, mask) {
                        Ok(Some(result)) => {
                            if result.1[0] > jd_start + ECLIPSE_WHEN_EPOCH_MARGIN_DAYS {
                                return Ok(result);
                            }
                        }
                        Ok(None) => {}
                        Err(err) => {
                            if !is_ephemeris_boundary(&err) {
                                return Err(err);
                            }
                        }
                    }
                }
            }
            Err(err) => {
                if !is_ephemeris_boundary(&err) {
                    return Err(err);
                }
            }
        }
    }

    let mut jd = jd_start;
    for _ in 0..ECLIPSE_WHEN_MAX_LUNATIONS {
        let jd_new_moon = {
            let mut state = |jd: f64| eph.sun_moon_state(jd);
            match find_next_new_moon(&mut state, jd) {
                Ok(jd_new_moon) => jd_new_moon,
                Err(err) => {
                    if is_ephemeris_boundary(&err) {
                        break;
                    }
                    return Err(err);
                }
            }
        };
        let mut result = match check_new_moon_for_eclipse(eph, jd_new_moon, mask) {
            Ok(result) => result,
            Err(err) => {
                if is_ephemeris_boundary(&err) {
                    break;
                }
                return Err(err);
            }
        };
        if let Some((_, times)) = &result {
            if times[0] <= jd_start + ECLIPSE_WHEN_EPOCH_MARGIN_DAYS {
                result = None;
            }
        }
        if let Some(found) = result {
            return Ok(found);
        }
        jd = jd_new_moon + ECLIPSE_WALK_ADVANCE_DAYS;
    }
    Err(Error::EclipseSearchExhausted {
        entry: "sol_eclipse_when_glob",
        detail: format!(
            "No matching solar eclipse found before the ephemeris boundary searching from JD {jd_start}"
        ),
    })
}

/// Public Native walk with type-mask admission before the first actual scalar
/// acquisition. It retains all lunation, bracket, contact and direction limits;
/// nonfinite clocks fail at the consumed position/time primitive. This variant
/// adds no eager finite-date guard to that acquisition order.
pub(crate) fn sol_eclipse_when_glob_source_via(
    eph: &impl SolarGlobalEphemeris,
    jd_start: f64,
    mask: i32,
    direction: EclipseSearchDirection,
) -> Result<(i32, [f64; 10])> {
    sol_glob_reject_impossible(mask)?;

    if direction == EclipseSearchDirection::Bidirectional {
        let jd_prev = {
            let mut state = |jd: f64| eph.sun_moon_state(jd);
            find_previous_new_moon(&mut state, jd_start)?
        };
        if jd_start - jd_prev <= ECLIPSE_PRE_PROBE_WINDOW_DAYS {
            if let Some(result) = check_new_moon_for_eclipse_consumed(eph, jd_prev, mask)? {
                if result.1[0] > jd_start + ECLIPSE_WHEN_EPOCH_MARGIN_DAYS {
                    return Ok(result);
                }
            }
        }
    }
    if direction == EclipseSearchDirection::Backward {
        let jd_next = {
            let mut state = |jd: f64| eph.sun_moon_state(jd);
            find_next_new_moon(&mut state, jd_start)?
        };
        if jd_next - jd_start <= ECLIPSE_PRE_PROBE_WINDOW_DAYS {
            if let Some(result) = check_new_moon_for_eclipse_consumed(eph, jd_next, mask)? {
                if result.1[0] < jd_start - ECLIPSE_WHEN_EPOCH_MARGIN_DAYS {
                    return Ok(result);
                }
            }
        }
        let mut jd = jd_start;
        for _ in 0..ECLIPSE_WHEN_MAX_LUNATIONS {
            let jd_prev = {
                let mut state = |jd: f64| eph.sun_moon_state(jd);
                match find_previous_new_moon(&mut state, jd) {
                    Ok(jd_prev) => jd_prev,
                    Err(err) => {
                        if is_ephemeris_boundary(&err) {
                            break;
                        }
                        return Err(err);
                    }
                }
            };
            match check_new_moon_for_eclipse_consumed(eph, jd_prev, mask) {
                Ok(Some(result)) => {
                    if result.1[0] < jd_start - ECLIPSE_WHEN_EPOCH_MARGIN_DAYS {
                        return Ok(result);
                    }
                }
                Ok(None) => {}
                Err(err) => {
                    if is_ephemeris_boundary(&err) {
                        break;
                    }
                    return Err(err);
                }
            }
            jd = jd_prev - 1.0;
        }
        return Err(Error::EclipseSearchExhausted {
            entry: "sol_eclipse_when_glob",
            detail: format!(
                "No matching solar eclipse found before the ephemeris boundary searching back from JD {jd_start}"
            ),
        });
    }

    if direction == EclipseSearchDirection::Forward {
        let previous = {
            let mut state = |jd: f64| eph.sun_moon_state(jd);
            match find_previous_new_moon(&mut state, jd_start) {
                Ok(jd) => Some(jd),
                Err(error) if is_ephemeris_boundary(&error) => None,
                Err(error) => return Err(error),
            }
        };
        if let Some(jd_prev) = previous {
            if jd_start - jd_prev <= ECLIPSE_PRE_PROBE_WINDOW_DAYS {
                if let Some(result) = check_new_moon_for_eclipse_consumed(eph, jd_prev, mask)? {
                    if result.1[0] > jd_start + ECLIPSE_WHEN_EPOCH_MARGIN_DAYS {
                        return Ok(result);
                    }
                }
            }
        }
    }

    let mut jd = jd_start;
    for _ in 0..ECLIPSE_WHEN_MAX_LUNATIONS {
        let jd_new_moon = {
            let mut state = |jd: f64| eph.sun_moon_state(jd);
            match find_next_new_moon(&mut state, jd) {
                Ok(jd_new_moon) => jd_new_moon,
                Err(err) => {
                    if is_ephemeris_boundary(&err) {
                        break;
                    }
                    return Err(err);
                }
            }
        };
        let mut result = match check_new_moon_for_eclipse_consumed(eph, jd_new_moon, mask) {
            Ok(result) => result,
            Err(err) => {
                if is_ephemeris_boundary(&err) {
                    break;
                }
                return Err(err);
            }
        };
        if let Some((_, times)) = &result {
            if times[0] <= jd_start + ECLIPSE_WHEN_EPOCH_MARGIN_DAYS {
                result = None;
            }
        }
        if let Some(found) = result {
            return Ok(found);
        }
        jd = jd_new_moon + ECLIPSE_WALK_ADVANCE_DAYS;
    }
    Err(Error::EclipseSearchExhausted {
        entry: "sol_eclipse_when_glob",
        detail: format!(
            "No matching solar eclipse found before the ephemeris boundary searching from JD {jd_start}"
        ),
    })
}

/// Selenocentric immersion depth at one instant, degrees:
/// the selenocentric angle between the solar and terrestrial limbs —
/// center separation less both apparent radii — minimized at maximum
/// eclipse. Negative while the Earth covers the Sun as seen from the
/// Moon. `moon_xyz_au` / `sun_xyz_au` are geocentric equatorial AU
/// vectors.
pub fn selenocentric_depth_deg(
    moon_xyz_au: (f64, f64, f64),
    sun_xyz_au: (f64, f64, f64),
) -> Result<f64> {
    for value in [
        moon_xyz_au.0,
        moon_xyz_au.1,
        moon_xyz_au.2,
        sun_xyz_au.0,
        sun_xyz_au.1,
        sun_xyz_au.2,
    ] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    let s = (
        sun_xyz_au.0 - moon_xyz_au.0,
        sun_xyz_au.1 - moon_xyz_au.1,
        sun_xyz_au.2 - moon_xyz_au.2,
    );
    let g = (-moon_xyz_au.0, -moon_xyz_au.1, -moon_xyz_au.2);
    let ds = (s.0.powi(2) + s.1.powi(2) + s.2.powi(2)).sqrt();
    let dm = (g.0.powi(2) + g.1.powi(2) + g.2.powi(2)).sqrt();
    if ds.is_nan() || ds <= 0.0 || dm.is_nan() || dm <= 0.0 {
        // Degenerate coincident vectors; a typed refusal instead of NaN words.
        return Err(Error::ZeroVectorNormalize);
    }
    let cosang = (s.0 * g.0 + s.1 * g.1 + s.2 * g.2) / (ds * dm);
    let ang = safe_acos(cosang).to_degrees();
    let rearth = (ECL_REARTH_AU / dm).asin().to_degrees();
    let rsun = (ECL_RSUN_AU / ds).asin().to_degrees();
    Ok(ang - (rearth + rsun))
}

/// Phase duration in minutes from two contact instants: `0.0` when either contact is unresolved.
pub fn eclipse_phase_duration_min(t_start_jd: f64, t_end_jd: f64) -> f64 {
    if t_start_jd == 0.0 || t_end_jd == 0.0 {
        0.0
    } else {
        (t_end_jd - t_start_jd) * 24.0 * 60.0
    }
}

/// Classify the eclipse (if any) at one Full Moon:
/// the coarse node pre-filter, selenocentric-depth maximum
/// refinement, canonical shadow classification and the type-mask test.
/// Shared by the pre-probe and the main walk.
pub(super) fn check_full_moon_for_eclipse<E: LunarGlobalEphemeris + ?Sized>(
    eph: &E,
    jd_full_moon: f64,
    mask: i32,
) -> Result<Option<(i32, [f64; 10])>> {
    let state = eph.sun_moon_state(jd_full_moon)?;
    let node_lon = eph.mean_node_lon_deg(jd_full_moon)?;
    let node_dist = moon_node_distance_deg(node_lon, state.moon_lon_deg)?;
    if node_dist >= LUNAR_NODE_PREFILTER_DEG {
        return Ok(None);
    }
    let jd_max = {
        let mut depth = |jd: f64| -> Result<f64> {
            let (moon_xyz, sun_xyz) = eph.geocentric_xyz_au(jd)?;
            selenocentric_depth_deg(moon_xyz, sun_xyz)
        };
        refine_lunar_eclipse_maximum(&mut depth, jd_full_moon)?
    };
    let (retc, _) = eph.lun_shadow(jd_max)?;
    if retc == 0 {
        return Ok(None);
    }
    let type_matches = (mask & ECL_TOTAL != 0 && retc & ECL_TOTAL != 0)
        || (mask & ECL_PARTIAL != 0 && retc & ECL_PARTIAL != 0)
        || (mask & ECL_PENUMBRAL != 0 && retc & ECL_PENUMBRAL != 0);
    if !type_matches {
        return Ok(None);
    }
    let mut shadow = |jd: f64| eph.lun_shadow(jd);
    let times = lun_eclipse_phase_times_via(jd_max, retc, &mut shadow)?;
    Ok(Some((retc, times)))
}

/// Reduce a flag word to its ephemeris-selection bits.
pub fn ecl_eph_flags(flags: i32) -> i32 {
    flags & (FLG_JPLEPH | FLG_SWIEPH | FLG_MOSEPH)
}

/// Local eclipse maximum from a separation provider: golden-section
/// minimization of the topocentric Sun-Moon separation over
/// `jd_approx ± search_range` (`60` iterations, `1e-8` day width
/// stop). Returns `(jd_maximum, min_separation_degrees)`.
pub fn calc_local_eclipse_max_time_via(
    separation_deg: &mut impl FnMut(f64) -> Result<f64>,
    jd_approx: f64,
    search_range_days: f64,
) -> Result<(f64, f64)> {
    if !jd_approx.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_approx });
    }
    if !search_range_days.is_finite() || search_range_days <= 0.0 {
        return Err(Error::InvalidSearchStep {
            value: search_range_days,
        });
    }
    let jd_max = golden_section_min(
        separation_deg,
        jd_approx - search_range_days,
        jd_approx + search_range_days,
        SOLAR_MAXIMUM_TOLERANCE_DAYS,
        GOLDEN_SECTION_ITERATIONS,
    )?;
    let min_separation = separation_deg(jd_max)?;
    if !min_separation.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: min_separation });
    }
    Ok((jd_max, min_separation))
}

/// Local eclipse maximum around a global maximum from a separation
/// provider. The bracket is resolved to the representable UT1-day spacing
/// around its center, with at most eighty golden-section contractions. This
/// is followed by eighth-order stationary-point refinement of the squared
/// separation at three fixed dyadic spacings. Squaring removes the narrow
/// square-root shape of a close angular passage; the minimum is unchanged for
/// a nonnegative separation. Search resolution is not a physical error bound.
/// Epochs and the positive half-window are UT1 days; separation is in degrees.
/// Returns the UT1 day of the minimum. Nonfinite bounds or a nonpositive window
/// are errors, and failures from any consumed provider sample propagate.
pub fn find_local_maximum_via(
    separation_deg: &mut impl FnMut(f64) -> Result<f64>,
    jd_center: f64,
    search_range_days: f64,
) -> Result<f64> {
    if !jd_center.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_center });
    }
    if !search_range_days.is_finite() || search_range_days <= 0.0 {
        return Err(Error::InvalidSearchStep {
            value: search_range_days,
        });
    }
    let mut minimum = golden_section_min(
        separation_deg,
        jd_center - search_range_days,
        jd_center + search_range_days,
        2.0 * (jd_center.next_up() - jd_center),
        2 * LOCAL_MAX_ITERATIONS,
    )?;
    // Differencing across a resolved interval avoids the nearly equal
    // objective values that limit a bracket-only solve near its minimum.
    // These spacings depend only on the declared search window, never on a
    // comparison residual, target epoch or returned circumstance attribute.
    for divisor in [128.0, 256.0, 512.0] {
        let step = search_range_days / divisor;
        let center = separation_deg(minimum)?.powi(2);
        let mut gradient = 0.0;
        let mut curvature = 0.0;
        for (index, (first, second)) in [
            (4.0 / 5.0, 8.0 / 5.0),
            (-1.0 / 5.0, -1.0 / 5.0),
            (4.0 / 105.0, 8.0 / 315.0),
            (-1.0 / 280.0, -1.0 / 560.0),
        ]
        .into_iter()
        .enumerate()
        {
            let offset = (index + 1) as f64 * step;
            let before = separation_deg(minimum - offset)?.powi(2);
            let after = separation_deg(minimum + offset)?.powi(2);
            gradient += first * (after - before);
            curvature += second * ((after - center) + (before - center));
        }
        if curvature > 0.0 {
            let correction = gradient / curvature * step;
            let next = minimum - correction;
            // A positive-curvature local refinement stays inside both its
            // resolved sampling interval and the original search bracket.
            if correction.abs() <= step
                && next > jd_center - search_range_days
                && next < jd_center + search_range_days
            {
                minimum = next;
            }
        }
    }
    Ok(minimum)
}

impl LocalEclipsePhases {
    /// All-zero (no local eclipse) circumstances.
    pub fn zeros() -> Self {
        Self { slots: [0.0; 10] }
    }

    /// Time of local maximum, JD UT (`slots[0]`; `0.0` when none).
    pub fn jd_max(&self) -> f64 {
        self.slots[0]
    }

    /// Eclipse magnitude at maximum (`slots[5]`).
    pub fn magnitude(&self) -> f64 {
        self.slots[5]
    }

    /// Borrow the ten slots in the release layout.
    pub fn as_array(&self) -> [f64; 10] {
        self.slots
    }
}

/// Local solar-eclipse circumstances at one global maximum: the `-1.0` deg Sun gate at the global maximum, the
/// `±3 h` / `30`-iteration golden search for the local minimum
/// separation, the Sun gate at the local maximum, the `(959.63 /
/// 3600)` / `(932.56 / 3600)` apparent-radius spellings, the
/// sum-of-radii overlap gate, the clamped magnitude plus the
/// unclamped diameter fraction plus the safe obscuration, the legacy
/// bisection contacts (`±2 h` outer, `±30 min` inner against the
/// absolute radius difference), and the horizon zeroing that skips
/// already-unresolved (`0.0`) contacts so no out-of-coverage epoch
/// is ever evaluated.
///
/// `sun_altaz_deg` answers `(true_altitude_deg, azimuth_deg)` and
/// `distances_au` answers `(sun_dist_au, moon_dist_au)`.
pub fn calculate_local_eclipse_phases_via(
    jd_max_global: f64,
    separation_deg: &mut impl FnMut(f64) -> Result<f64>,
    sun_altaz_deg: &mut impl FnMut(f64) -> Result<(f64, f64)>,
    distances_au: &mut impl FnMut(f64) -> Result<(f64, f64)>,
) -> Result<LocalEclipsePhases> {
    if !jd_max_global.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: jd_max_global,
        });
    }
    let (sun_alt_global, _) = sun_altaz_deg(jd_max_global)?;
    if sun_alt_global < LOCAL_SUN_ALT_GATE_DEG {
        return Ok(LocalEclipsePhases::zeros());
    }
    let jd_local_max = golden_section_min(
        separation_deg,
        jd_max_global - LOCAL_MAX_SEARCH_RANGE_DAYS,
        jd_max_global + LOCAL_MAX_SEARCH_RANGE_DAYS,
        LOCAL_PHASES_TOLERANCE_DAYS,
        LOCAL_PHASES_GOLDEN_ITERATIONS,
    )?;
    let min_separation = separation_deg(jd_local_max)?;
    if !min_separation.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: min_separation });
    }
    let (sun_alt_max, sun_az_max) = sun_altaz_deg(jd_local_max)?;
    if sun_alt_max < LOCAL_SUN_ALT_GATE_DEG {
        return Ok(LocalEclipsePhases::zeros());
    }
    let (sun_dist_au, moon_dist_au) = distances_au(jd_local_max)?;
    for value in [sun_dist_au, moon_dist_au] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    let sun_angular_radius = (SUN_SEMIDIAMETER_ARCSEC / 3600.0) / sun_dist_au;
    let moon_angular_radius =
        (MOON_SEMIDIAMETER_ARCSEC / 3600.0) * (MOON_SEMIDIAMETER_REF_AU / moon_dist_au);
    let sun_diameter = 2.0 * sun_angular_radius;
    let moon_diameter = 2.0 * moon_angular_radius;
    let sum_radii = sun_angular_radius + moon_angular_radius;
    if min_separation > sum_radii {
        return Ok(LocalEclipsePhases::zeros());
    }
    let overlap = sum_radii - min_separation;
    let magnitude = if sun_diameter < MINIMUM_SEPARATION_FOR_LENS_DEG {
        0.0
    } else {
        (overlap / sun_diameter).clamp(0.0, 1.0 + moon_diameter / sun_diameter)
    };
    let fraction_covered = if sun_diameter < MINIMUM_SEPARATION_FOR_LENS_DEG {
        0.0
    } else {
        // Release: the [6] slot keeps the unclamped fraction
        // while [5] carries the clamped magnitude.
        overlap / sun_diameter
    };
    let obscuration =
        calculate_obscuration(sun_angular_radius, moon_angular_radius, min_separation)?;
    let mut jd_first = find_local_contact_time(
        separation_deg,
        jd_local_max - LOCAL_CONTACT_SEARCH_DAYS,
        jd_local_max,
        sum_radii,
        false,
    )?;
    let mut jd_fourth = find_local_contact_time(
        separation_deg,
        jd_local_max,
        jd_local_max + LOCAL_CONTACT_SEARCH_DAYS,
        sum_radii,
        true,
    )?;
    // Horizon zeroing skips unresolved (`0.0`) contacts.
    if jd_first != 0.0 && sun_altaz_deg(jd_first)?.0 < LOCAL_SUN_ALT_GATE_DEG {
        jd_first = 0.0;
    }
    if jd_fourth != 0.0 && sun_altaz_deg(jd_fourth)?.0 < LOCAL_SUN_ALT_GATE_DEG {
        jd_fourth = 0.0;
    }
    let mut jd_second = 0.0;
    let mut jd_third = 0.0;
    let diff_radii = (moon_angular_radius - sun_angular_radius).abs();
    if min_separation < diff_radii {
        jd_second = find_local_contact_time(
            separation_deg,
            jd_local_max - LOCAL_INNER_CONTACT_SEARCH_DAYS,
            jd_local_max,
            diff_radii,
            false,
        )?;
        jd_third = find_local_contact_time(
            separation_deg,
            jd_local_max,
            jd_local_max + LOCAL_INNER_CONTACT_SEARCH_DAYS,
            diff_radii,
            true,
        )?;
        if jd_second != 0.0 && sun_altaz_deg(jd_second)?.0 < LOCAL_SUN_ALT_GATE_DEG {
            jd_second = 0.0;
        }
        if jd_third != 0.0 && sun_altaz_deg(jd_third)?.0 < LOCAL_SUN_ALT_GATE_DEG {
            jd_third = 0.0;
        }
    }
    Ok(LocalEclipsePhases {
        slots: [
            jd_local_max,
            jd_first,
            jd_second,
            jd_third,
            jd_fourth,
            magnitude,
            fraction_covered,
            obscuration,
            sun_az_max,
            sun_alt_max,
        ],
    })
}

/// Local eclipse type flags: `0` when nothing is covered, else `ECL_VISIBLE` plus
/// `ECL_CENTRAL` with `ECL_TOTAL` / `ECL_ANNULAR_TOTAL` /
/// `ECL_ANNULAR` for a central eclipse at the location, or
/// `ECL_PARTIAL` otherwise.
pub fn determine_local_eclipse_type(magnitude: f64, moon_sun_ratio: f64, is_central: bool) -> i32 {
    if magnitude <= 0.0 {
        return 0;
    }
    let mut flags = ECL_VISIBLE;
    if is_central {
        flags |= ECL_CENTRAL;
        if moon_sun_ratio >= 1.0 {
            flags |= ECL_TOTAL;
        } else if moon_sun_ratio > 0.99 {
            flags |= ECL_ANNULAR_TOTAL;
        } else {
            flags |= ECL_ANNULAR;
        }
    } else {
        flags |= ECL_PARTIAL;
    }
    flags
}

/// Instantaneous solar-eclipse magnitude at a site: `0.0` when the Sun's true altitude is below the
/// `-1.0` deg gate or the discs do not overlap, else the
/// diameter-cover fraction clamped to `[0, 1 + moon/sun]`.
pub fn sol_eclipse_magnitude_at_loc_via(
    sun_alt_true_deg: f64,
    separation_deg: f64,
    sun_dist_au: f64,
    moon_dist_au: f64,
) -> Result<f64> {
    for value in [sun_alt_true_deg, separation_deg, sun_dist_au, moon_dist_au] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    if sun_alt_true_deg < LOCAL_SUN_ALT_GATE_DEG {
        return Ok(0.0);
    }
    let sun_angular_radius = (SUN_SEMIDIAMETER_ARCSEC / 3600.0) / sun_dist_au;
    let moon_angular_radius =
        (MOON_SEMIDIAMETER_ARCSEC / 3600.0) * (MOON_SEMIDIAMETER_REF_AU / moon_dist_au);
    let sum_radii = sun_angular_radius + moon_angular_radius;
    if separation_deg >= sum_radii {
        return Ok(0.0);
    }
    let sun_diameter = 2.0 * sun_angular_radius;
    let moon_diameter = 2.0 * moon_angular_radius;
    let overlap = sum_radii - separation_deg;
    Ok((overlap / sun_diameter).clamp(0.0, 1.0 + moon_diameter / sun_diameter))
}

/// Next solar eclipse visible from a site: the global-lunation walk with the
/// first-iteration `_WHEN_LOC_EPOCH_WINDOW` offer, the served
/// topocentric-maximum search, the full-disc overlap gate, the
/// annular/total/partial phase rule, the outer (`±2 h`, `±4 h`
/// retry) and inner (`±30 min`) contact solves, the local-contact
/// epoch gate (plus the re-entrancy arm so `jd = tret[0];
/// when_loc(jd)` advances), the apparent-altitude visibility bits,
/// the sunrise/sunset clamping with the horizon-crossed maximum
/// move, and the `how`/`where` circumstance assembly with
/// `dcore[0]` in `attr[3]`.
///
/// Returns `(eclipse_type, tret10, attr20)`: `tret` holds maximum /
/// first / second / third / fourth contact plus sunrise/sunset in
/// `[5]`/`[6]`; `attr` is the 20-slot `how` array. Geographic
/// coordinates live provider-side, so the exhaustion detail omits longitude
/// and latitude. Any forward global-candidate error becomes search exhaustion;
/// any backward candidate error ends the current collection sweep, retaining
/// candidates already found. Local-maximum coverage errors become exhaustion;
/// other failures in local geometry, contacts and circumstances propagate.
pub fn sol_eclipse_when_loc_via(
    eph: &impl SolarLocalEphemeris,
    tjdut: f64,
    backwards: bool,
) -> Result<(i32, [f64; 10], [f64; 20])> {
    if !tjdut.is_finite() {
        return Err(Error::NonFiniteArgument { value: tjdut });
    }
    let mut jd = tjdut;
    let mut first_iteration = true;
    for _ in 0..SOL_LOC_MAX_ECLIPSES {
        // Next/previous global candidate. Forward:
        // the current lunation is offered even when its geocentric
        // maximum just slipped behind the epoch (first iteration
        // only); the local-contact gate below keeps or skips it.
        // Backward: collect forward from a year back and keep the
        // most recent eclipse before the epoch window.
        let jd_max_global = if backwards {
            let earlier_jd = jd - SOL_LOC_BACKWARD_COLLECT_DAYS;
            let horizon = jd + WHEN_LOC_EPOCH_WINDOW_DAYS;
            let mut found: Option<(i32, [f64; 10])> = None;
            let mut temp_jd = earlier_jd;
            while temp_jd < horizon {
                match eph.next_global_forward(temp_jd) {
                    Ok((global_type, global_times)) => {
                        if global_times[0] < horizon {
                            found = Some((global_type, global_times));
                        }
                        temp_jd = global_times[0] + ECLIPSE_WALK_ADVANCE_DAYS;
                    }
                    Err(_) => break,
                }
            }
            match found {
                Some((_, global_times)) => global_times[0],
                None => {
                    jd -= SOL_LOC_BACKWARD_COLLECT_DAYS;
                    continue;
                }
            }
        } else {
            let anchor = if first_iteration {
                jd - WHEN_LOC_EPOCH_WINDOW_DAYS
            } else {
                jd
            };
            first_iteration = false;
            match eph.next_global_forward(anchor) {
                Ok((_, global_times)) => global_times[0],
                Err(_) => {
                    return Err(Error::EclipseSearchExhausted {
                        entry: "sol_eclipse_when_loc",
                        detail: format!(
                            "No solar eclipse visible from the observer within 50 years of JD {tjdut}"
                        ),
                    });
                }
            }
        };

        let mut separation = |jd: f64| eph.separation_deg(jd);
        let jd_local_max = match find_local_maximum_via(
            &mut separation,
            jd_max_global,
            LOCAL_MAX_SEARCH_RANGE_DAYS,
        ) {
            Ok(jd_local_max) => jd_local_max,
            Err(err) => {
                if is_ephemeris_boundary(&err) {
                    return Err(Error::EclipseSearchExhausted {
                        entry: "sol_eclipse_when_loc",
                        detail: format!(
                            "No solar eclipse visible from the observer within 50 years of JD {tjdut}"
                        ),
                    });
                }
                return Err(err);
            }
        };
        let (rsun, rmoon, min_separation) = eph.radii_sep_deg(jd_local_max, false)?;
        if min_separation > rsun + rmoon {
            jd = if backwards {
                jd_max_global - 1.0
            } else {
                jd_max_global + ECLIPSE_WALK_ADVANCE_DAYS
            };
            continue;
        }
        // Local phase at maximum.
        let phase = if min_separation < rsun - rmoon {
            ECL_ANNULAR
        } else if min_separation < (rsun - rmoon).abs() {
            ECL_TOTAL
        } else {
            ECL_PARTIAL
        };
        let two_hours = 2.0 / 24.0;
        let mut f_outer = |jd_c: f64| -> Result<f64> {
            let (rs_c, rm_c, sep_c) = eph.radii_sep_deg(jd_c, false)?;
            Ok(sep_c - (rs_c + rm_c))
        };
        let mut f_inner = |jd_c: f64| -> Result<f64> {
            let (rs_c, rm_c, sep_c) = eph.radii_sep_deg(jd_c, true)?;
            Ok((rs_c - rm_c).abs() - sep_c)
        };
        let mut jd_first =
            bisect_sign_change(&mut f_outer, jd_local_max - two_hours, jd_local_max)?;
        if jd_first == 0.0 {
            jd_first =
                bisect_sign_change(&mut f_outer, jd_local_max - 2.0 * two_hours, jd_local_max)?;
        }
        let mut jd_fourth =
            bisect_sign_change(&mut f_outer, jd_local_max, jd_local_max + two_hours)?;
        if jd_fourth == 0.0 {
            jd_fourth =
                bisect_sign_change(&mut f_outer, jd_local_max, jd_local_max + 2.0 * two_hours)?;
        }
        // Epoch gate on the LOCAL contacts: an
        // eclipse whose local phases are still unfolding at the
        // epoch is kept; a wholly past (forward) or wholly future
        // (backward) one is skipped. The second arm stops the
        // re-entrant `jd = tret[0]; when_loc(jd)` idiom from
        // re-finding the same eclipse.
        let first_contact = if jd_first != 0.0 {
            jd_first
        } else {
            jd_local_max
        };
        let last_contact = if jd_fourth != 0.0 {
            jd_fourth
        } else {
            jd_local_max
        };
        if (!backwards
            && (last_contact <= tjdut + ECLIPSE_WHEN_EPOCH_MARGIN_DAYS
                || jd_local_max <= tjdut + ECLIPSE_WHEN_EPOCH_MARGIN_DAYS))
            || (backwards
                && (first_contact >= tjdut - ECLIPSE_WHEN_EPOCH_MARGIN_DAYS
                    || jd_local_max >= tjdut - ECLIPSE_WHEN_EPOCH_MARGIN_DAYS))
        {
            jd = if backwards {
                jd_max_global - 1.0
            } else {
                jd_max_global + ECLIPSE_WALK_ADVANCE_DAYS
            };
            continue;
        }
        let mut jd_second = 0.0;
        let mut jd_third = 0.0;
        if min_separation < (rsun - rmoon).abs() {
            let half_hour = LOCAL_INNER_CONTACT_SEARCH_DAYS;
            jd_second = bisect_sign_change(&mut f_inner, jd_local_max - half_hour, jd_local_max)?;
            jd_third = bisect_sign_change(&mut f_inner, jd_local_max, jd_local_max + half_hour)?;
        }
        // Visibility bits: a phase counts when the
        // Sun's APPARENT altitude is positive there; `tret` keeps the
        // geometric times regardless.
        let mut ecl_type = phase;
        for (tc, bit) in [
            (jd_local_max, ECL_MAX_VISIBLE),
            (jd_first, ECL_1ST_VISIBLE),
            (jd_second, ECL_2ND_VISIBLE),
            (jd_third, ECL_3RD_VISIBLE),
            (jd_fourth, ECL_4TH_VISIBLE),
        ] {
            if tc == 0.0 {
                continue;
            }
            if eph.sun_altaz_deg(tc)?.2 > 0.0 {
                ecl_type |= ECL_VISIBLE | bit;
            }
        }
        if ecl_type & ECL_VISIBLE == 0 {
            jd = if backwards {
                jd_max_global - 1.0
            } else {
                jd_max_global + ECLIPSE_WALK_ADVANCE_DAYS
            };
            continue;
        }
        // Sunrise/sunset bounding the visible part: a circumpolar Sun (`-2`) leaves the geometric
        // times untouched; a wholly nocturnal eclipse advances the
        // walk; a horizon crossing inside the contacts is reported
        // and, when the geometric maximum itself is down, becomes
        // the reported maximum.
        let mut jd_reported_max = jd_local_max;
        let mut jd_sunrise = 0.0;
        let mut jd_sunset = 0.0;
        let mut clamped = false;
        let crossings = eph.sun_horizon_crossings(first_contact - 0.001)?;
        let circumpolar = crossings.rise_rc == -2 || crossings.set_rc == -2;
        if !circumpolar {
            let tjdr = crossings.rise_jd;
            let tjds = crossings.set_jd;
            if tjds < jd_first || (tjds > tjdr && tjdr > jd_fourth) {
                jd = if backwards {
                    jd_max_global - 1.0
                } else {
                    jd_max_global + ECLIPSE_WALK_ADVANCE_DAYS
                };
                continue;
            }
            if jd_first < tjdr && tjdr < jd_fourth {
                jd_sunrise = tjdr;
                if ecl_type & ECL_MAX_VISIBLE == 0 {
                    jd_reported_max = tjdr;
                    clamped = true;
                }
            }
            if jd_first < tjds && tjds < jd_fourth {
                jd_sunset = tjds;
                if ecl_type & ECL_MAX_VISIBLE == 0 {
                    jd_reported_max = tjds;
                    clamped = true;
                }
            }
        }
        let (retc_how, mut attr) = eph.how_core(jd_reported_max)?;
        if clamped {
            ecl_type &= !(ECL_TOTAL | ECL_ANNULAR | ECL_PARTIAL);
            ecl_type |= retc_how & (ECL_TOTAL | ECL_ANNULAR | ECL_PARTIAL);
        }
        let (rc_where, dcore) = eph.where_core(jd_reported_max)?;
        ecl_type |= rc_where & ECL_NONCENTRAL;
        attr[3] = dcore[0];
        let tret = [
            jd_reported_max,
            jd_first,
            jd_second,
            jd_third,
            jd_fourth,
            jd_sunrise,
            jd_sunset,
            0.0,
            0.0,
            0.0,
        ];
        return Ok((ecl_type, tret, attr));
    }
    Err(Error::EclipseSearchExhausted {
        entry: "sol_eclipse_when_loc",
        detail: format!("No solar eclipse visible from the observer within 50 years of JD {tjdut}"),
    })
}

/// Physical equatorial disc radius of an occulted body in AU: star names (`None`) and unknown ids
/// are point sources with radius `0.0`; the Sun and Moon use the
/// [`ECL_RSUN_AU`] / [`ECL_RMOON_AU`] identities; planets use the release
/// `_PLANET_RADIUS_KM` table and the remaining finite bodies (Chiron,
/// Pholus, Ceres, Pallas, Juno, Vesta, Earth) the `_EXTRA_BODY_RADIUS_KM`
/// table, both over [`ECL_AU_KM`].
pub fn occ_body_radius_au(body_id: Option<i32>) -> f64 {
    let id = match body_id {
        None => return 0.0,
        Some(id) => id,
    };
    if id == SUN {
        return ECL_RSUN_AU;
    }
    if id == MOON {
        return ECL_RMOON_AU;
    }
    let radius_km = match id {
        MERCURY => Some(2439.4),
        VENUS => Some(6051.8),
        MARS => Some(3389.5),
        JUPITER => Some(69911.0),
        SATURN => Some(58232.0),
        URANUS => Some(25362.0),
        NEPTUNE => Some(24622.0),
        PLUTO => Some(1188.3),
        CHIRON => Some(83.0),
        PHOLUS => Some(95.0),
        CERES => Some(469.7),
        PALLAS => Some(256.5),
        JUNO => Some(123.298),
        VESTA => Some(261.385),
        EARTH => Some(6371.0084),
        _ => None,
    };
    match radius_km {
        Some(km) => km / ECL_AU_KM,
        None => 0.0,
    }
}

/// Dynamic angular disc radius of a planet in degrees from its geocentric
/// distance in AU: `asin(min(1, radius_km / dist_km))` over the
/// planet-radius table only; the extra finite bodies in
/// [`occ_body_radius_au`] are point sources here. Unknown IDs and non-positive
/// distances yield [`OCC_POINT_SOURCE_RADIUS_DEG`]. For a recognized planet,
/// positive infinite distance yields zero radius; a NaN distance reaches
/// `f64::min(1.0)` and produces a right-angle radius rather than propagating
/// NaN. This raw geometric helper does not validate a physical distance.
pub fn calc_planet_angular_radius_deg(planet_id: i32, dist_au: f64) -> f64 {
    let radius_km = match planet_id {
        MERCURY => Some(2439.4),
        VENUS => Some(6051.8),
        MARS => Some(3389.5),
        JUPITER => Some(69911.0),
        SATURN => Some(58232.0),
        URANUS => Some(25362.0),
        NEPTUNE => Some(24622.0),
        PLUTO => Some(1188.3),
        _ => None,
    };
    let radius_km = match radius_km {
        None => return OCC_POINT_SOURCE_RADIUS_DEG,
        Some(km) => km,
    };
    let dist_km = dist_au * ECL_AU_KM;
    if dist_km <= 0.0 {
        return OCC_POINT_SOURCE_RADIUS_DEG;
    }
    (radius_km / dist_km).min(1.0).asin().to_degrees()
}

/// Answer a negative body id as the Sun.
pub fn fold_pseudo_body_to_sun(body_id: i32) -> i32 {
    if body_id < 0 { SUN } else { body_id }
}

/// Point-source contact coincidence.
/// Infallible tuple — times travel.
pub fn collapse_point_source_contacts(jd_second: f64, jd_third: f64) -> (f64, f64) {
    (jd_second, jd_third)
}

/// Advance `t` to the Moon-body conjunction in ecliptic longitude
/// first-order Newton iteration on the longitude gap with
/// [`OCC_MOON_RATE_DEG_DAY`], converging to [`OCC_CONJ_TOL_DEG`] within
/// [`OCC_CONJUNCTION_MAX_STEPS`] steps. The first gap is `(blon - mlon)
/// % 360` with only the backward `-360` shift; later gaps fold past
/// 180 deg, exactly as in the release.
///
/// Returns `Ok(None)` when the FIRST body evaluation steps off the
/// active ephemeris. Every other provider failure —
/// including failures during the Newton iteration itself — propagates
/// unchanged, matching the release `try` placement. A star past
/// [`OCC_STAR_LAT_LIMIT_DEG`] is [`Error::InvalidEclipseSearch`] (the
/// release base-`Error`; the serving binding formats the star name into
/// the diagnostic detail).
///
/// Returns `(t_conj, body_ecl_lat, moon_ecl_lat)` at the converged
/// instant. `direction` is negative for a backward search.
pub fn seek_moon_conjunction(
    body_lonlat_deg: &mut impl FnMut(f64) -> Result<(f64, f64)>,
    moon_lonlat_deg: &mut impl FnMut(f64) -> Result<(f64, f64)>,
    t_start: f64,
    direction: i32,
    is_star: bool,
) -> Result<Option<(f64, f64, f64)>> {
    if !t_start.is_finite() {
        return Err(Error::NonFiniteArgument { value: t_start });
    }
    let (mut blon, mut blat) = match body_lonlat_deg(t_start) {
        Ok(pair) => pair,
        Err(err) => {
            if is_ephemeris_boundary(&err) {
                return Ok(None);
            }
            return Err(err);
        }
    };
    if is_star && star_never_occulted(blat) {
        return Err(Error::InvalidEclipseSearch {
            entry: "lun_occult",
            detail: format!("occultation never occurs: star has ecl. lat. {blat:.1}"),
        });
    }
    let mut mpos = moon_lonlat_deg(t_start)?;
    let mut gap = (blon - mpos.0).rem_euclid(360.0);
    if direction < 0 {
        gap -= 360.0;
    }
    let mut t = t_start;
    let mut steps = 0usize;
    while gap.abs() > OCC_CONJ_TOL_DEG && steps < OCC_CONJUNCTION_MAX_STEPS {
        t += gap / OCC_MOON_RATE_DEG_DAY;
        (blon, blat) = body_lonlat_deg(t)?;
        mpos = moon_lonlat_deg(t)?;
        gap = (blon - mpos.0).rem_euclid(360.0);
        if gap > 180.0 {
            gap -= 360.0;
        }
        steps += 1;
    }
    Ok(Some((t, blat, mpos.1)))
}

/// Wrap a raw body-minus-Moon right-ascension difference into
/// `(-180, 180]` degrees.
pub fn wrap_ra_gap_deg(raw_diff_deg: f64) -> f64 {
    let mut d = raw_diff_deg.rem_euclid(360.0);
    if d > 180.0 {
        d -= 360.0;
    }
    d
}

/// `(azimuth_deg, true_altitude_deg, apparent_altitude_deg)` for a
/// topocentric ecliptic place.
/// Nutation follows the release `azalt` path (the full IAU 2006/2000A
/// series, here the native [`crate::erfa::nut06a`]), not the stored
/// LEB series the calc path evaluates: the release computes these two
/// legs from different nutation evaluations, so each leg matches its
/// own counterpart. The frame epoch rides the bundled default-timescale
/// TT; the sidereal leg
/// keeps the caller-model Delta-T, exactly like [`crate::events::azalt_deg`].
/// Shared by the lunar-local serving and the planetary-occultation
/// local serving so the `azalt`-equivalent spelling exists once.
#[allow(clippy::too_many_arguments)]
pub(super) fn altaz_of_topo_ecliptic_impl(
    topo_ecl: (f64, f64, f64),
    jd_ut: f64,
    lon_deg_east: f64,
    lat_deg: f64,
    alt_m: f64,
    deltat_at: &dyn Fn(f64) -> Result<f64>,
    offset_lo_deg: f64,
    offset_hi_deg: f64,
) -> Result<(f64, f64, f64)> {
    use crate::time::{SiderealDeltaT, local_apparent_sidereal_time_hours};
    if !jd_ut.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_ut });
    }
    for v in [topo_ecl.0, topo_ecl.1, topo_ecl.2] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    let dt_days_first = deltat_at(jd_ut)?;
    if !dt_days_first.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: dt_days_first,
        });
    }
    // Frame epoch for nutation/obliquity: the bundled default-timescale TT
    // not the caller-model TT.
    let jd_tt_frame = crate::time::default_timescale_tt_from_ut1(jd_ut)?;
    let (dpsi_rad, deps_rad) = crate::erfa::nut06a(jd_tt_frame)?;
    let eps_true_deg =
        crate::positions::vondrak_mean_obliquity_deg(jd_tt_frame)? + deps_rad.to_degrees();
    let (ra_deg, dec_deg, _) = cotrans(topo_ecl.0, topo_ecl.1, topo_ecl.2, -eps_true_deg);
    let dt_days = deltat_at(jd_ut)?;
    if !dt_days.is_finite() {
        return Err(Error::NonFiniteArgument { value: dt_days });
    }
    let sdt = SiderealDeltaT::try_new(dt_days, offset_lo_deg, offset_hi_deg)?;
    let last_hours = local_apparent_sidereal_time_hours(
        jd_ut,
        lon_deg_east,
        eps_true_deg,
        dpsi_rad.to_degrees(),
        &sdt,
    )?;
    let ha_deg = (last_hours * 15.0 - ra_deg).rem_euclid(360.0);
    let ha_rad = ha_deg.to_radians();
    let dec_rad = dec_deg.to_radians();
    let lat_rad = lat_deg.to_radians();
    let sin_alt = lat_rad.sin() * dec_rad.sin() + lat_rad.cos() * dec_rad.cos() * ha_rad.cos();
    let true_alt_deg = sin_alt.clamp(-1.0, 1.0).asin().to_degrees();
    let y = ha_rad.sin() * dec_rad.cos();
    let x = ha_rad.cos() * dec_rad.cos() * lat_rad.sin() - dec_rad.sin() * lat_rad.cos();
    let azimuth_deg = y.atan2(x).to_degrees().rem_euclid(360.0);
    let pressure_mbar = crate::refraction::estimate_pressure_mbar(alt_m);
    let (apparent_alt_deg, _) = crate::refraction::refrac_extended_true_to_app(
        true_alt_deg,
        alt_m,
        pressure_mbar,
        10.0,
        RISE_SET_DIP_LAPSE_RATE_K_M,
    );
    if !apparent_alt_deg.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: apparent_alt_deg,
        });
    }
    Ok((azimuth_deg, true_alt_deg, apparent_alt_deg))
}

/// Validate an eclipse path sampling interval and convert it to days
/// following the release `_eclipse_sampling_step_days` gate order, with
/// Typed errors: non-finite window bounds are [`Error::EphemerisRange`]
/// with the same semantics as `validate_jd_range`; step refusals are
/// [`Error::InvalidEclipseSearch`] carrying the release message in
/// `detail` under the caller's `func_name` entry).
///
/// The resolution test is exact: doubling a finite double is exact and
/// the doubles spacing is a power of two, so requiring the step to be
/// strictly more than half the spacing at the window's largest-magnitude
/// instant guarantees the `jd += step_days` accumulator advances at every
/// visited sample. Half the spacing exactly is refused (round-half-even
/// advances only from odd mantissas).
pub fn eclipse_sampling_step_days(
    jd_start: f64,
    jd_end: f64,
    step_minutes: f64,
    func_name: &'static str,
) -> Result<f64> {
    if !jd_start.is_finite() || !jd_end.is_finite() {
        return Err(Error::EphemerisRange {
            requested_jd: if !jd_start.is_finite() {
                jd_start
            } else {
                jd_end
            },
            start_jd: jd_start,
            end_jd: jd_end,
            body: None,
        });
    }
    // Checked before the sign test: `nan <= 0.0` is false, so a NaN step
    // would slip through the positivity test below.
    if !step_minutes.is_finite() {
        return Err(Error::InvalidEclipseSearch {
            entry: func_name,
            detail: format!(
                "{func_name}: step_minutes must be a finite number of minutes, got {step_minutes:?}"
            ),
        });
    }
    // `<= 0.0` refuses `0.0` and `-0.0` as well: neither advances the
    // accumulator.
    if step_minutes <= 0.0 {
        return Err(Error::InvalidEclipseSearch {
            entry: func_name,
            detail: format!(
                "{func_name}: step_minutes must be greater than zero, got {step_minutes:?}"
            ),
        });
    }
    let step_days = step_minutes / MINUTES_PER_DAY;
    let pivot = jd_start.abs().max(jd_end.abs());
    let spacing = pivot.next_up() - pivot;
    if 2.0 * step_days <= spacing {
        let min_minutes = (spacing / 2.0) * MINUTES_PER_DAY;
        return Err(Error::InvalidEclipseSearch {
            entry: func_name,
            detail: format!(
                "{func_name}: step_minutes={step_minutes:?} cannot advance the \
                 sampling instant near Julian Day {pivot:.1}, where consecutive \
                 double-precision values are {spacing:.3e} days apart; use a step \
                 larger than {min_minutes:.3e} minutes"
            ),
        });
    }
    Ok(step_days)
}

/// Parse a geographic `(lon, lat, alt)` triple with the same length rule
/// as the release `geopos` handling in `sol_eclipse_when_loc` /
/// `_sol_eclipse_how_details_impl`: at least 3 elements, `(lon, lat, alt)`
/// order, extras ignored.
/// `message` is one of [`GEOPOS_SEQUENCE_MESSAGE`] /
/// [`GEOPOS_ELEMENTS_MESSAGE`] per caller; non-finite components travel
/// untouched exactly as the release `float` conversion passes them, and
/// range validation stays with the consumer (`validate_geopos`).
pub fn parse_geopos(
    values: &[f64],
    entry: &'static str,
    message: &'static str,
) -> Result<(f64, f64, f64)> {
    if values.len() < 3 {
        return Err(Error::InvalidEclipseSearch {
            entry,
            detail: message.to_string(),
        });
    }
    Ok((values[0], values[1], values[2]))
}

/// Normalize a path longitude to `[-180, 180]`.
pub fn normalize_path_longitude_deg(lon_deg: f64) -> f64 {
    (lon_deg + 180.0).rem_euclid(360.0) - 180.0
}

/// Geographic points along the central line of a solar eclipse in the
/// same walk order as the release `calc_eclipse_central_line`: validated
/// step, `while jd <= jd_end` accumulation, `ECL_CENTRAL` gate,
/// `gamma < 1.5` surface gate, nonzero-coordinates gate; partial-only
/// eclipses contribute no points.
/// Returns `(times_jd_ut, latitudes_deg, longitudes_deg)`; provider
/// failures propagate.
///
/// Sampled points carry the underlying where-composition precision of about
/// 0.001° (roughly 100 m on the ground). The shadow geometry follows Meeus,
/// *Astronomical Algorithms*, chapter 54; published-path comparisons use the
/// Espenak & Meeus *Five Millennium Canon of Solar Eclipses*; background
/// theory is chapter 11 of the *Explanatory Supplement to the Astronomical
/// Almanac*. The 2024-04-08 total eclipse is the canonical worked example
/// for this walk.
pub fn eclipse_central_line_via(
    eph: &impl EclipsePathEphemeris,
    jd_start: f64,
    jd_end: f64,
    step_minutes: f64,
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>)> {
    let step_days =
        eclipse_sampling_step_days(jd_start, jd_end, step_minutes, "calc_eclipse_central_line")?;
    let mut times = Vec::new();
    let mut latitudes = Vec::new();
    let mut longitudes = Vec::new();
    let mut jd = jd_start;
    while jd <= jd_end {
        let (retflag, lon, lat) = eph.path_where(jd)?;
        if retflag & ECL_CENTRAL != 0 {
            let (x, y) = eph.path_besselian_xy(jd)?;
            let gamma = (x * x + y * y).sqrt();
            if gamma < PATH_CENTRAL_GAMMA_LIMIT
                && (lon.abs() > PATH_CENTRAL_COORD_EPS_DEG
                    || lat.abs() > PATH_CENTRAL_COORD_EPS_DEG)
            {
                times.push(jd);
                latitudes.push(lat);
                longitudes.push(lon);
            }
        }
        jd += step_days;
    }
    Ok((times, latitudes, longitudes))
}

/// One northern/southern limit point in `(latitude, longitude)` degrees,
/// following the release limit-sampler operation order: `|l2|` umbral
/// radius with the `0.001` floor, the `gamma < 1 + f + umbra` touch gate,
/// the `y ± |l2|` edge offset — north adds, south subtracts — the surface
/// gate, the `z`-factor arm, clamped latitude, the polar `cos_lat` floor,
/// the `-mu + atan2` longitude, the exact
/// `tan(phi_geod) = tan(phi_geoc) / (1 - f)^2` conversion and floored
/// normalization). `None` means the sample is omitted (central shadow
/// misses Earth here); non-finite `x`/`y`/`l2` words also omit via the
/// comparison gates. A NaN declination propagates to a NaN latitude here
pub fn eclipse_limit_point_deg(
    x_er: f64,
    y_er: f64,
    d_deg: f64,
    mu_deg: f64,
    l2_er: f64,
    south: bool,
) -> Option<(f64, f64)> {
    let umbra_radius = l2_er.abs();
    let gamma = (x_er * x_er + y_er * y_er).sqrt();
    if !(gamma < 1.0 + PATH_LIMIT_FLATTENING + umbra_radius
        && umbra_radius > PATH_LIMIT_UMBRA_FLOOR_ER)
    {
        return None;
    }
    let y_edge = if south {
        y_er - umbra_radius
    } else {
        y_er + umbra_radius
    };
    let gamma_edge = (x_er * x_er + y_edge * y_edge).sqrt();
    // Surface gate via a binding rather than `!(a < b)`: NaN compares
    // false, so a non-finite edge is omitted exactly as the release
    // comparison gate omits it.
    let on_surface = gamma_edge < 1.0 + PATH_LIMIT_FLATTENING;
    if !on_surface {
        return None;
    }
    let z_factor = if gamma_edge > PATH_LIMIT_SURFACE_GAMMA {
        0.0
    } else {
        (0.0f64.max(1.0 - gamma_edge * gamma_edge)).sqrt()
    };
    let d_rad = d_deg.to_radians();
    let (sin_d, cos_d) = (d_rad.sin(), d_rad.cos());
    let sin_lat = (y_edge * cos_d + z_factor * sin_d).clamp(-1.0, 1.0);
    let lat = sin_lat.asin().to_degrees();
    let cos_lat = lat.to_radians().cos();
    let lon_offset = if cos_lat > PATH_LIMIT_COS_LAT_FLOOR {
        (x_er.atan2(z_factor * cos_d - y_edge * sin_d)).to_degrees()
    } else {
        0.0
    };
    let lon = -mu_deg + lon_offset;
    let lat_geodetic = (lat.to_radians().tan() / (1.0 - PATH_LIMIT_FLATTENING).powi(2))
        .atan()
        .to_degrees();
    Some((lat_geodetic, normalize_path_longitude_deg(lon)))
}

/// Shared limit walk with the same loop shape as the release
/// northern/southern loops: validated step, `while jd <= jd_end`
/// accumulation, samples with no surface point omitted. `south` selects
/// the edge; `func_name` owns the step errors.
pub(super) fn eclipse_limit_walk_via(
    eph: &impl EclipseLimitEphemeris,
    south: bool,
    jd_start: f64,
    jd_end: f64,
    step_minutes: f64,
    func_name: &'static str,
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>)> {
    let step_days = eclipse_sampling_step_days(jd_start, jd_end, step_minutes, func_name)?;
    let mut times = Vec::new();
    let mut latitudes = Vec::new();
    let mut longitudes = Vec::new();
    let mut jd = jd_start;
    while jd <= jd_end {
        let (x, y, d, mu, l2) = eph.limit_besselian(jd)?;
        if let Some((lat, lon)) = eclipse_limit_point_deg(x, y, d, mu, l2, south) {
            times.push(jd);
            latitudes.push(lat);
            longitudes.push(lon);
        }
        jd += step_days;
    }
    Ok((times, latitudes, longitudes))
}

/// Points along the northern limit of a central eclipse path (same
/// construction as the release `calc_eclipse_northern_limit`, over
/// [`eclipse_limit_point_deg`]).
///
/// Accuracy degrades near the grazing ends of a central path, where the
/// shadow meets the surface tangentially. References: Meeus, *Astronomical
/// Algorithms*, chapter 54; *Explanatory Supplement to the Astronomical
/// Almanac*, chapter 11.
pub fn eclipse_northern_limit_via(
    eph: &impl EclipseLimitEphemeris,
    jd_start: f64,
    jd_end: f64,
    step_minutes: f64,
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>)> {
    eclipse_limit_walk_via(
        eph,
        false,
        jd_start,
        jd_end,
        step_minutes,
        "calc_eclipse_northern_limit",
    )
}

/// Points along the southern limit of a central eclipse path (same
/// construction as the release `calc_eclipse_southern_limit`, over
/// [`eclipse_limit_point_deg`]).
///
/// Accuracy degrades near the grazing ends of a central path, where the
/// shadow meets the surface tangentially. References: Meeus, *Astronomical
/// Algorithms*, chapter 54; *Explanatory Supplement to the Astronomical
/// Almanac*, chapter 11.
pub fn eclipse_southern_limit_via(
    eph: &impl EclipseLimitEphemeris,
    jd_start: f64,
    jd_end: f64,
    step_minutes: f64,
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>)> {
    eclipse_limit_walk_via(
        eph,
        true,
        jd_start,
        jd_end,
        step_minutes,
        "calc_eclipse_southern_limit",
    )
}

/// Classify the outer error variant for a leaf eclipse computation.
/// Direct LEB/SPK corruption and coverage errors propagate. Other variants,
/// including errors wrapping those causes, follow the existing zero/skip
/// degradation policy; this classifier does not inspect `primary_error()`.
pub(super) fn eclipse_leaf_reraises(err: &Error) -> bool {
    matches!(
        err,
        Error::EphemerisRange { .. }
            | Error::LebOutOfRange { .. }
            | Error::SpkOutOfRange { .. }
            | Error::LebCorrupt { .. }
            | Error::SpkCorrupt { .. }
    )
}

/// Width of the totality/annularity path in km in the same order as the
/// release `_calc_eclipse_path_width_impl`: explicit-site vs central-line
/// location, `ECL_CENTRAL` gate on the where arm only, Sun-altitude gate,
/// disc-overlap gate, the umbral/antumbral cone geometry, the `sin_alt`
/// floor and the 1500 km cap). `lat_deg`/`lon_deg` are both-or-neither
/// (a half-provided site falls back to the central line, as the release
/// `lat is not None and lon is not None` test does).
///
/// Simplified-model caveat: the cone radius is projected onto the surface
/// through `1/sin` of the Sun altitude, measuring the shadow in the local
/// vertical plane rather than perpendicular to the ground track, so near
/// greatest eclipse (Sun high) widths read a few percent above the published
/// canon values (about 5–8% high on the banked comparisons: roughly 124 vs
/// 115 km for 2017-08-21 and 207 vs 198 km for 2024-04-08). The gap grows
/// toward sunrise/sunset where grazing geometry dominates; the model assumes
/// a spherical Earth with the WGS84 equatorial radius. References: Meeus,
/// *Astronomical Algorithms*, chapter 54; Espenak & Meeus *Five Millennium
/// Canon of Solar Eclipses*. The 2024-04-08 total eclipse (Dallas site) is
/// the canonical worked example for this computation.
///
/// # Errors
///
/// Non-finite request epochs, central-line acquisition errors and checked
/// sample-domain errors propagate. A failed `path_width_sample` propagates only
/// direct range/corruption variants selected by the leaf classifier; other
/// failures, including wrapped causes, currently degrade to `Ok(0.0)`.
pub fn eclipse_path_width_via(
    eph: &impl EclipsePathWidthEphemeris,
    jd_ut: f64,
    lat_deg: Option<f64>,
    lon_deg: Option<f64>,
) -> Result<f64> {
    if !jd_ut.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_ut });
    }
    let (calc_lon, calc_lat) = match (lat_deg, lon_deg) {
        (Some(lat), Some(lon)) => (lon, lat),
        _ => {
            let (retflag, lon, lat) = eph.path_width_where(jd_ut)?;
            if retflag & ECL_CENTRAL == 0 {
                return Ok(0.0);
            }
            (lon, lat)
        }
    };
    let sample = match eph.path_width_sample(jd_ut, calc_lon, calc_lat) {
        Ok(sample) => sample,
        // The existing leaf policy inspects only the outer error variant.
        // Other sample failures, including wrapped causes and horizontal
        // reductions, retain their zero-width degradation here.
        Err(e) if eclipse_leaf_reraises(&e) => return Err(e),
        Err(_) => return Ok(0.0),
    };
    if sample.sun_alt_true_deg <= 0.0 {
        return Ok(0.0);
    }
    for value in [
        sample.sun_dist_au,
        sample.moon_dist_au,
        sample.sun_alt_true_deg,
        sample.separation_deg,
    ] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    let (sun_angular_radius, moon_angular_radius) =
        eclipse_disc_angular_radii_deg(sample.sun_dist_au, sample.moon_dist_au)?;
    let moon_sun_ratio = moon_angular_radius / sun_angular_radius;
    let diff_radii = (moon_angular_radius - sun_angular_radius).abs();
    if sample.separation_deg > diff_radii {
        return Ok(0.0);
    }
    let sun_dist_km = sample.sun_dist_au * ECL_AU_KM;
    let moon_dist_km = sample.moon_dist_au * ECL_AU_KM;
    let alpha = ((PATH_WIDTH_SUN_RADIUS_KM - PATH_WIDTH_MOON_RADIUS_KM) / sun_dist_km).atan();
    let umbra_radius_km = if moon_sun_ratio >= 1.0 {
        let umbra = PATH_WIDTH_MOON_RADIUS_KM - moon_dist_km * alpha.tan();
        if umbra <= 0.0 {
            return Ok(0.0);
        }
        umbra
    } else {
        let umbra = moon_dist_km * alpha.tan() - PATH_WIDTH_MOON_RADIUS_KM;
        if umbra < 0.0 { umbra.abs() } else { umbra }
    };
    let mut sin_alt = sample.sun_alt_true_deg.to_radians().sin();
    if sin_alt < PATH_WIDTH_SIN_ALT_FLOOR {
        sin_alt = PATH_WIDTH_SIN_ALT_FLOOR;
    }
    let path_width_km = (2.0 * umbra_radius_km / sin_alt).min(PATH_WIDTH_MAX_KM);
    Ok(path_width_km.max(0.0))
}

/// Central-line maximum duration of totality/annularity in minutes,
/// with the release `calc_solar_eclipse_duration` composition:
/// non-central eclipses and missing C2/C3 yield `0.0`; otherwise
/// `(C3 - C2)` in minutes at the shadow-center point — not the global
/// U1/U4 span.
///
/// Central durations match published local-maximum durations to about ±0.05
/// minutes. References: Meeus, *Astronomical Algorithms*, chapter 54;
/// Espenak & Meeus *Five Millennium Canon of Solar Eclipses*; *Explanatory
/// Supplement to the Astronomical Almanac* (2013), chapter 11. The
/// 2024-04-08 total eclipse is the canonical worked example.
pub fn solar_eclipse_duration_min_via(
    eph: &impl SolarDurationEphemeris,
    jd_max: f64,
) -> Result<f64> {
    if !jd_max.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_max });
    }
    let (retflag, center_lon, center_lat) = eph.duration_where(jd_max)?;
    if retflag & ECL_CENTRAL == 0 {
        return Ok(0.0);
    }
    let (jd_c2, jd_c3) = eph.duration_local_contacts(jd_max, center_lat, center_lon)?;
    if jd_c2 <= 0.0 || jd_c3 <= 0.0 {
        return Ok(0.0);
    }
    Ok((jd_c3 - jd_c2) * MINUTES_PER_DAY)
}

/// Precise maximum-eclipse time dispatcher with the release
/// `sol_eclipse_max_time` composition: `None` finds the global gamma
/// minimum, `Some((lat, lon))` the local separation minimum via
/// [`calc_local_eclipse_max_time_via`]. Both-or-neither lat/lon is enforced by the `Option` type; a non-finite epoch is refused and a
/// non-positive search window is [`Error::InvalidSearchStep`], mirroring
/// the local core's gates.
///
/// Greatest-eclipse positions from this composition are good to about 0.02°.
pub fn eclipse_max_time_via(
    gamma: &mut impl FnMut(f64) -> Result<f64>,
    separation_deg: &mut impl FnMut(f64) -> Result<f64>,
    jd_approx: f64,
    site_lat_lon_deg: Option<(f64, f64)>,
    search_range_days: f64,
) -> Result<(f64, f64)> {
    if !jd_approx.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_approx });
    }
    if !search_range_days.is_finite() || search_range_days <= 0.0 {
        return Err(Error::InvalidSearchStep {
            value: search_range_days,
        });
    }
    match site_lat_lon_deg {
        None => {
            let jd_max = golden_section_min(
                gamma,
                jd_approx - search_range_days,
                jd_approx + search_range_days,
                SOLAR_MAXIMUM_TOLERANCE_DAYS,
                GOLDEN_SECTION_ITERATIONS,
            )?;
            let gamma_at_max = gamma(jd_max)?;
            if !gamma_at_max.is_finite() {
                return Err(Error::NonFiniteSearchValue { x: gamma_at_max });
            }
            Ok((jd_max, gamma_at_max))
        }
        Some(_) => calc_local_eclipse_max_time_via(separation_deg, jd_approx, search_range_days),
    }
}

/// Kernel coverage used to clamp occultation batch scans with the same
/// end-selection rule as the release `_kernel_jd_bounds`: minimum segment
/// start, maximum segment end — the maximum, not the minimum, so split
/// kernels keep their true upper bound. Empty input yields `None` (the
/// release `min` raises); segment order is irrelevant.
pub fn kernel_jd_bounds(segments_start_end_jd: &[(f64, f64)]) -> Option<(f64, f64)> {
    let mut iter = segments_start_end_jd.iter();
    let (mut earliest, mut latest) = *iter.next()?;
    for &(start, end) in iter {
        if start < earliest {
            earliest = start;
        }
        if end > latest {
            latest = end;
        }
    }
    Some((earliest, latest))
}

/// Position angle of the Moon relative to the Sun in degrees with the
/// release `_calc_position_angle` formula: `atan2` of the RA/Dec pair,
/// folded to `[0, 360)`. Inputs are equatorial degrees; non-finite words
/// are refused.
pub fn position_angle_deg(
    sun_ra_deg: f64,
    sun_dec_deg: f64,
    moon_ra_deg: f64,
    moon_dec_deg: f64,
) -> Result<f64> {
    for value in [sun_ra_deg, sun_dec_deg, moon_ra_deg, moon_dec_deg] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    let sun_ra = sun_ra_deg.to_radians();
    let sun_dec = sun_dec_deg.to_radians();
    let delta_ra = moon_ra_deg.to_radians() - sun_ra;
    let y = delta_ra.sin();
    let x = sun_dec.cos() * moon_dec_deg.to_radians().tan() - sun_dec.sin() * delta_ra.cos();
    let mut pa = y.atan2(x).to_degrees();
    if pa < 0.0 {
        pa += 360.0;
    }
    Ok(pa)
}

/// Maximum magnitude at the local maximum with the release inline
/// composition spelling: `max(0.0, min(overlap / sun_diameter, 1.0 +
/// moon_diameter / sun_diameter))`.
pub fn details_max_magnitude(sun_diameter_deg: f64, moon_diameter_deg: f64, sep_deg: f64) -> f64 {
    let overlap = (sun_diameter_deg / 2.0 + moon_diameter_deg / 2.0) - sep_deg;
    0.0f64.max((overlap / sun_diameter_deg).min(1.0 + moon_diameter_deg / sun_diameter_deg))
}

/// Comprehensive local solar-eclipse circumstances in the release
/// `_sol_eclipse_how_details_impl` composition order: initial-type
/// assembly, early return when no eclipse is visible, the 50-iteration
/// local maximum, the maximum magnitude/obscuration/type update, the
/// `±2.5 h` outer and `±2.5/4 h` inner contacts each isolated so one
/// non-range provider failure keeps the remaining components (LEB
/// coverage misses propagate — release re-raise), and the durations).
///
/// `(lon_deg, lat_deg, alt_m)` is the observer site. The site travels
/// inside the provider closures (as the release `geopos` does); it is
/// listed so the serving binding threads it to acquisition, and the
/// `geopos` length rule is the caller's ([`parse_geopos`]).
///
/// The shadow model follows Meeus, *Astronomical Algorithms*, chapter 54;
/// the 2024-04-08 total eclipse is the canonical worked example for local
/// solar circumstances.
pub fn eclipse_how_details_via(
    eph: &impl EclipseHowDetailsEphemeris,
    jd_ut: f64,
    _lon_deg: f64,
    _lat_deg: f64,
    _alt_m: f64,
) -> Result<EclipseHowDetails> {
    if !jd_ut.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_ut });
    }
    let (eclipse_type, attr) = eph.details_initial(jd_ut)?;
    let mut out = EclipseHowDetails {
        eclipse_type,
        is_visible: eclipse_type & ECL_VISIBLE != 0,
        is_total: eclipse_type & ECL_TOTAL != 0,
        is_annular: eclipse_type & ECL_ANNULAR != 0,
        is_partial: eclipse_type & ECL_PARTIAL != 0,
        jd_c1: 0.0,
        jd_c2: 0.0,
        jd_max: 0.0,
        jd_c3: 0.0,
        jd_c4: 0.0,
        magnitude: attr[0],
        max_magnitude: 0.0,
        obscuration: attr[2],
        max_obscuration: 0.0,
        obscuration_percent: attr[2] * 100.0,
        max_obscuration_percent: 0.0,
        ratio: attr[1],
        shadow_width_km: attr[3],
        position_angle_c1: 0.0,
        position_angle_c2: 0.0,
        position_angle_c3: 0.0,
        position_angle_c4: 0.0,
        sun_alt_c1: 0.0,
        sun_az_c1: 0.0,
        sun_alt_c2: 0.0,
        sun_az_c2: 0.0,
        sun_alt_max: attr[5],
        sun_az_max: attr[4],
        sun_alt_c3: 0.0,
        sun_az_c3: 0.0,
        sun_alt_c4: 0.0,
        sun_az_c4: 0.0,
        duration_partial_minutes: 0.0,
        duration_total_minutes: 0.0,
        sun_angular_radius: 0.0,
        moon_angular_radius: 0.0,
        separation: attr[7],
    };
    if eclipse_type == 0 {
        return Ok(out);
    }
    let (sun_r, moon_r) = eph.details_disc_radii_deg(jd_ut)?;
    out.sun_angular_radius = sun_r;
    out.moon_angular_radius = moon_r;
    let diff_radii = (sun_r - moon_r).abs();
    let mut separation = |jd: f64| eph.details_separation_deg(jd);
    let jd_local_max = golden_section_min(
        &mut separation,
        jd_ut - DETAILS_LOCAL_MAX_RANGE_DAYS,
        jd_ut + DETAILS_LOCAL_MAX_RANGE_DAYS,
        DETAILS_LOCAL_MAX_TOL_DAYS,
        DETAILS_LOCAL_MAX_ITERATIONS,
    )?;
    out.jd_max = jd_local_max;
    let max_sep = eph.details_separation_deg(jd_local_max)?;
    let (sun_alt_max, sun_az_max) = eph.details_sun_alt_az_deg(jd_local_max)?;
    out.sun_alt_max = sun_alt_max;
    out.sun_az_max = sun_az_max;
    let (sun_r_max, moon_r_max) = eph.details_disc_radii_deg(jd_local_max)?;
    let sun_diameter = 2.0 * sun_r_max;
    let moon_diameter = 2.0 * moon_r_max;
    out.max_magnitude = details_max_magnitude(sun_diameter, moon_diameter, max_sep);
    out.ratio = moon_diameter / sun_diameter;
    let obscuration = super::geometry::details_max_obscuration(sun_r_max, moon_r_max, max_sep);
    let max_obscuration = obscuration.max(0.0);
    out.max_obscuration = max_obscuration;
    out.max_obscuration_percent = max_obscuration * 100.0;
    let is_central = max_sep < diff_radii;
    if is_central {
        out.is_partial = false;
        if moon_r_max >= sun_r_max {
            out.is_total = true;
            out.is_annular = false;
        } else {
            out.is_total = false;
            out.is_annular = true;
        }
    } else {
        out.is_partial = true;
        out.is_total = false;
        out.is_annular = false;
    }
    let contact_range = DETAILS_CONTACT_RANGE_DAYS;
    if let Some(block) = details_contact_block(
        eph,
        &mut separation,
        jd_local_max - contact_range,
        jd_local_max,
        sun_r_max + moon_r_max,
        false,
    )? {
        out.jd_c1 = block.jd;
        if let Some((alt, az)) = block.alt_az {
            out.sun_alt_c1 = alt;
            out.sun_az_c1 = az;
        }
        if let Some(pa) = block.position_angle_deg {
            out.position_angle_c1 = pa;
        }
    }
    if let Some(block) = details_contact_block(
        eph,
        &mut separation,
        jd_local_max,
        jd_local_max + contact_range,
        sun_r_max + moon_r_max,
        true,
    )? {
        out.jd_c4 = block.jd;
        if let Some((alt, az)) = block.alt_az {
            out.sun_alt_c4 = alt;
            out.sun_az_c4 = az;
        }
        if let Some(pa) = block.position_angle_deg {
            out.position_angle_c4 = pa;
        }
    }
    let diff_radii_max = (sun_r_max - moon_r_max).abs();
    if is_central && max_sep < diff_radii_max {
        if let Some(block) = details_contact_block(
            eph,
            &mut separation,
            jd_local_max - contact_range / 4.0,
            jd_local_max,
            diff_radii_max,
            false,
        )? {
            out.jd_c2 = block.jd;
            if let Some((alt, az)) = block.alt_az {
                out.sun_alt_c2 = alt;
                out.sun_az_c2 = az;
            }
            if let Some(pa) = block.position_angle_deg {
                out.position_angle_c2 = pa;
            }
        }
        if let Some(block) = details_contact_block(
            eph,
            &mut separation,
            jd_local_max,
            jd_local_max + contact_range / 4.0,
            diff_radii_max,
            true,
        )? {
            out.jd_c3 = block.jd;
            if let Some((alt, az)) = block.alt_az {
                out.sun_alt_c3 = alt;
                out.sun_az_c3 = az;
            }
            if let Some(pa) = block.position_angle_deg {
                out.position_angle_c3 = pa;
            }
        }
    }
    if out.jd_c1 > 0.0 && out.jd_c4 > 0.0 {
        out.duration_partial_minutes = (out.jd_c4 - out.jd_c1) * MINUTES_PER_DAY;
    }
    if out.jd_c2 > 0.0 && out.jd_c3 > 0.0 {
        out.duration_total_minutes = (out.jd_c3 - out.jd_c2) * MINUTES_PER_DAY;
    }
    Ok(out)
}

/// Rotate one ecliptic or equatorial direction, preserving its distance.
///
/// Longitude, latitude and obliquity `eps_deg` are degrees. Negative obliquity
/// maps ecliptic longitude/latitude to right ascension/declination; positive
/// obliquity maps back. The rise/set altitude engine uses negative of-date
/// true obliquity for its ecliptic-to-equatorial rotation. The Cartesian
/// rotation and non-finite behavior are described in [`cotrans_radians`].
/// Longitude is reduced modulo 360 degrees;
/// zero is positive, while floating-point rounding can produce 360 degrees
/// at the upper endpoint. Distance passes through unchanged. The latitude
/// sine is clamped to `[-1, 1]`; NaN propagates through the clamp and `asin`.
/// Longitude also retains IEEE NaN propagation.
pub fn cotrans(lon_deg: f64, lat_deg: f64, dist: f64, eps_deg: f64) -> (f64, f64, f64) {
    let lon_rad = lon_deg.to_radians();
    let lat_rad = lat_deg.to_radians();
    let eps_rad = (-eps_deg).to_radians();
    cotrans_radians(lon_rad, lat_rad, eps_rad, dist)
}

/// Rotate a spherical direction using already evaluated radians.
///
/// Longitude, latitude and rotation angle `eps_rad` are radians. For a
/// degree-based obliquity convention, the rotation angle is `radians(-eps)`.
/// Returns longitude reduced modulo 360 degrees and latitude in degrees,
/// followed by the unchanged distance in the caller's units. Floating-point
/// rounding can produce a longitude of 360 degrees at the upper endpoint.
/// With `c = cos(latitude)` and rotation angle `e = eps_rad`, the rotated
/// components are `x = c cos(longitude)`,
/// `y = c sin(longitude) cos(e) - sin(latitude) sin(e)` and
/// `z = sin(latitude) cos(e) + c sin(e) sin(longitude)`.
/// Longitude is recovered with `atan2(y, x)` and latitude with `asin(z)`.
/// Full Cartesian components avoid division by the input latitude cosine,
/// including beyond the poles.
/// No conversion, finite-input check or provider acquisition is performed.
/// The latitude-sine clamp preserves NaN, so an unordered sine yields NaN
/// latitude; longitude also retains IEEE NaN propagation. A zero longitude has
/// the positive sign; latitude retains its computed sign. This is a closed-form rotation.
pub fn cotrans_radians(lon_rad: f64, lat_rad: f64, eps_rad: f64, dist: f64) -> (f64, f64, f64) {
    let cos_eps = eps_rad.cos();
    let sin_eps = eps_rad.sin();
    let cos_lat = lat_rad.cos();
    let sin_lat = lat_rad.sin();
    let cos_lon = lon_rad.cos();
    let sin_lon = lon_rad.sin();
    let sin_new_lat = sin_lat * cos_eps + cos_lat * sin_eps * sin_lon;
    // Bound finite rounding excursions without replacing a NaN sine.
    let new_lat_rad = sin_new_lat.clamp(-1.0, 1.0).asin();
    // Full-3D components (no `1 / cos_lat` division), so the longitude
    // retains the correct quadrant beyond the conventional latitude interval.
    let y = cos_lat * sin_lon * cos_eps - sin_lat * sin_eps;
    let x = cos_lat * cos_lon;
    let new_lon_rad = y.atan2(x);
    let new_lon = new_lon_rad.to_degrees().rem_euclid(360.0);
    (
        // Canonicalize only longitude zero, preserving nonzero remainders.
        if new_lon == 0.0 { 0.0 } else { new_lon },
        new_lat_rad.to_degrees(),
        dist,
    )
}

/// Composition family for a fictitious id, or `None` for a recognised id
/// with no served analytical model (Nibiru 49, the unnamed 59-69 tail).
/// Waldemath 58 is served through the release anchor derivation below
/// (Delta-T plus the sealed LEB nutation pair, both bound at the
/// serving).
pub(super) fn fictitious_family(body_id: i32) -> Option<FictFamily> {
    use crate::constants::{
        CUPIDO, HARRINGTON, ISIS, NEPTUNE_ADAMS, NEPTUNE_LEVERRIER, PLUTO_LOWELL, PLUTO_PICKERING,
        POSEIDON, PROSERPINA, VULCAN, WALDEMATH, WHITE_MOON,
    };
    if (CUPIDO..=POSEIDON).contains(&body_id) {
        Some(FictFamily::Uranian)
    } else if body_id == ISIS {
        Some(FictFamily::Transpluto)
    } else if body_id == HARRINGTON
        || body_id == NEPTUNE_LEVERRIER
        || body_id == NEPTUNE_ADAMS
        || body_id == PLUTO_LOWELL
        || body_id == PLUTO_PICKERING
        || body_id == VULCAN
        || body_id == PROSERPINA
    {
        Some(FictFamily::PredictedHelio)
    } else if body_id == WHITE_MOON {
        Some(FictFamily::WhiteMoon)
    } else if body_id == WALDEMATH {
        Some(FictFamily::Waldemath)
    } else {
        None
    }
}

/// True for the fictitious ids served on this path: the 18 model-backed
/// ids 40-48 and 50-58. Model-less recognised ids (Nibiru 49, the 59-69
/// tail) are refused at construction with the release probe's
/// illegal-body mapping.
pub fn is_fictitious_event_body(body_id: i32) -> bool {
    fictitious_family(body_id).is_some()
}

/// First-order annual aberration shift (same `xg += evel * rg / C`
/// association as the release `_geo_of_date` / Transpluto / Selena arms):
/// the shift of the geometric direction by the observer's barycentric
/// velocity over c. Always applied on the served apparent path.
pub(super) fn first_order_aberration_xyz(geo: [f64; 3], vel: [f64; 3], rg: f64) -> [f64; 3] {
    use crate::positions::C_LIGHT_AU_DAY;
    [
        geo[0] + vel[0] * rg / C_LIGHT_AU_DAY,
        geo[1] + vel[1] * rg / C_LIGHT_AU_DAY,
        geo[2] + vel[2] * rg / C_LIGHT_AU_DAY,
    ]
}

/// Unit-direction annual aberration.
pub(super) fn annual_aberration_cartesian(dir: [f64; 3], vel: [f64; 3]) -> [f64; 3] {
    use crate::positions::C_LIGHT_AU_DAY;
    let (rx, ry, rz) = (dir[0], dir[1], dir[2]);
    let r_mag = (rx * rx + ry * ry + rz * rz).sqrt();
    if r_mag < 1e-30 {
        return [0.0, 0.0, 0.0];
    }
    let (rx, ry, rz) = (rx / r_mag, ry / r_mag, rz / r_mag);
    let (vx, vy, vz) = (
        vel[0] / C_LIGHT_AU_DAY,
        vel[1] / C_LIGHT_AU_DAY,
        vel[2] / C_LIGHT_AU_DAY,
    );
    let r_dot_v = rx * vx + ry * vy + rz * vz;
    [vx - r_dot_v * rx, vy - r_dot_v * ry, vz - r_dot_v * rz]
}

/// Bradley annual aberration with renormalization (same form as the
/// release `astrometry.apply_aberration_to_position`, incl. both
/// `1e-30` gates): the Uranian geocentric branch correction.
pub(super) fn bradley_aberration_xyz(pos: [f64; 3], vel: [f64; 3]) -> [f64; 3] {
    let (x, y, z) = (pos[0], pos[1], pos[2]);
    let distance = (x * x + y * y + z * z).sqrt();
    if distance < 1e-30 {
        return pos;
    }
    let (ux, uy, uz) = (x / distance, y / distance, z / distance);
    let (dx, dy, dz) = {
        let d = annual_aberration_cartesian([ux, uy, uz], vel);
        (d[0], d[1], d[2])
    };
    let (ax, ay, az) = (ux + dx, uy + dy, uz + dz);
    let a_mag = (ax * ax + ay * ay + az * az).sqrt();
    if a_mag < 1e-30 {
        return pos;
    }
    [
        ax / a_mag * distance,
        ay / a_mag * distance,
        az / a_mag * distance,
    ]
}

/// Ecliptic Cartesian to spherical `(lon [0, 360), lat, dist)` (same
/// geocentric-branch shaping as the release: `% 360` longitude,
/// clamped `asin` latitude with the `r == 0` zero-latitude arm).
pub(super) fn ecl_xyz_to_spherical(v: [f64; 3]) -> (f64, f64, f64) {
    let rg = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    let lon = v[1].atan2(v[0]).to_degrees().rem_euclid(360.0);
    let lat = if rg > 0.0 {
        (v[2] / rg).clamp(-1.0, 1.0).asin().to_degrees()
    } else {
        0.0
    };
    (lon, lat, rg)
}

/// Mean-to-true ecliptic longitude.
pub(super) fn add_of_date_nutation_lon(
    serving: &crate::positions::SealedServing<'_, impl crate::positions::StateProvider>,
    lon_deg: f64,
    jd_tt: f64,
) -> Result<f64> {
    let (dpsi_rad, _) = (serving.nutation_provider())(jd_tt)?;
    Ok((lon_deg + dpsi_rad.to_degrees()).rem_euclid(360.0))
}

/// Magnitude-preserving ecliptic-J2000→of-date vector rotation (same
/// `_apply_hypothetical_topocentric._rotate_to_of_date` behavior as the
/// release: the zero vector maps to zeros; otherwise the direction
/// travels through the J2000→of-date precession with its magnitude kept).
pub(super) fn rotate_vector_to_of_date(v: [f64; 3], jd_tt: f64) -> Result<[f64; 3]> {
    let mag = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if mag == 0.0 {
        return Ok([0.0, 0.0, 0.0]);
    }
    let lon_v = v[1].atan2(v[0]).to_degrees();
    let lat_v = (v[2] / mag).clamp(-1.0, 1.0).asin().to_degrees();
    let (lon_d, lat_d) =
        crate::hypothetical::precess_ecliptic(lon_v, lat_v, crate::hypothetical::J2000_JD, jd_tt)?;
    let lr = lon_d.to_radians();
    let br = lat_d.to_radians();
    let cl = br.cos();
    Ok([mag * cl * lr.cos(), mag * cl * lr.sin(), mag * br.sin()])
}

// Sealed saros-series + solar max-time servings

/// Sealed global solar maximum time near `jd_approx`: the release
/// `sol_eclipse_max_time` global leg (golden-section minimum of gamma
/// over `jd_approx ± search_range_days`, 60 iterations, `1e-8` day
/// width stop). Returns `(jd_maximum, gamma_at_maximum)` with gamma in
/// Earth radii, re-read at the solved instant exactly as the release
/// re-evaluates `_calc_gamma(jd_max)`.
pub fn sol_global_max_time_sealed<P: crate::positions::StateProvider>(
    serving: &SealedSolarGlobalServing<'_, P>,
    jd_approx: f64,
    search_range_days: f64,
) -> Result<(f64, f64)> {
    let mut gamma_at = |jd: f64| serving.gamma_er(jd);
    let jd_max = refine_solar_eclipse_maximum(&mut gamma_at, jd_approx, search_range_days)?;
    Ok((jd_max, gamma_at(jd_max)?))
}

/// Classification with separate consumed where roles. true requests near-miss
/// magnitude; false requests central/umbral bits only. The same source algebra,
/// finite-domain checks and hybrid ordering apply; unused callbacks never run.
pub(crate) fn classify_solar_eclipse_consumed(
    state: &SunMoonState,
    gamma: f64,
    l1_er: f64,
    l2_er: f64,
    where_gate: &mut impl FnMut(bool) -> Result<(i32, f64)>,
    hybrid_probe: &mut impl FnMut() -> Result<bool>,
) -> Result<SolarEclipseClassification> {
    for value in [gamma, l1_er, l2_er] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    classify_solar_eclipse_stages(state, gamma, l1_er, l2_er, where_gate, hybrid_probe)
}

/// Classification over admitted words, with a branch-specific where request.
/// `true` consumes near-miss magnitude; `false` consumes centrality. Hybrid
/// acquisition remains behind the same core-shadow and centrality guards.
fn classify_solar_eclipse_stages(
    state: &SunMoonState,
    gamma: f64,
    l1_er: f64,
    l2_er: f64,
    where_gate: &mut impl FnMut(bool) -> Result<(i32, f64)>,
    hybrid_probe: &mut impl FnMut() -> Result<bool>,
) -> Result<SolarEclipseClassification> {
    // Angular radii in degrees.
    let sun_angular_radius = (959.63 / 3600.0) / state.sun_dist_au;
    let moon_angular_radius = (932.56 / 3600.0) * (MOON_SEMIDIAMETER_REF_AU / state.moon_dist_au);
    let moon_sun_ratio = if sun_angular_radius < 1e-10 {
        1.0
    } else {
        moon_angular_radius / sun_angular_radius
    };
    let gamma_limit_partial = 1.0 + l1_er;

    let separation = angular_separation_ecliptic(
        state.sun_lon_deg,
        state.sun_lat_deg,
        state.moon_lon_deg,
        state.moon_lat_deg,
    );
    let sum_radii = sun_angular_radius + moon_angular_radius;
    let mut magnitude = if separation >= sum_radii {
        0.0
    } else {
        ((sum_radii - separation) / (2.0 * sun_angular_radius)).clamp(0.0, ECLIPSE_MAGNITUDE_MAX)
    };

    if gamma.abs() > gamma_limit_partial {
        return Ok(SolarEclipseClassification {
            eclipse_type: 0,
            magnitude: 0.0,
            gamma,
            moon_sun_ratio,
        });
    }
    if is_near_miss_eclipse_with_limit(gamma, gamma_limit_partial) {
        let (_, where_magnitude) = where_gate(true)?;
        if where_magnitude <= 0.0 {
            return Ok(SolarEclipseClassification {
                eclipse_type: 0,
                magnitude: 0.0,
                gamma,
                moon_sun_ratio,
            });
        }
        let eclipse_type = ECL_PARTIAL | ECL_NONCENTRAL;
        magnitude =
            calculate_eclipse_magnitude_with_limit(gamma, moon_sun_ratio, gamma_limit_partial)?;
        return Ok(SolarEclipseClassification {
            eclipse_type,
            magnitude,
            gamma,
            moon_sun_ratio,
        });
    }

    let mut eclipse_type = ECL_PARTIAL | ECL_NONCENTRAL;
    if gamma.abs() <= 1.0 + l2_er.abs() && l2_er.abs() > 0.0 {
        let (where_retc, _) = where_gate(false)?;
        let central_bit = if where_retc & ECL_CENTRAL != 0 {
            ECL_CENTRAL
        } else if (where_retc & ECL_NONCENTRAL) != 0 && (where_retc & ECL_PARTIAL) == 0 {
            ECL_NONCENTRAL
        } else {
            0
        };
        if central_bit != 0 {
            if hybrid_probe()? {
                eclipse_type = ECL_ANNULAR_TOTAL | central_bit;
            } else if l2_er < 0.0 {
                eclipse_type = ECL_TOTAL | central_bit;
            } else {
                eclipse_type = ECL_ANNULAR | central_bit;
            }
        }
    }
    Ok(SolarEclipseClassification {
        eclipse_type,
        magnitude,
        gamma,
        moon_sun_ratio,
    })
}

/// Native candidate classification consumes only its requested where branch.
/// The prefilter, maximum, hybrid/contact arithmetic and type gate remain exact.
fn check_new_moon_for_eclipse_consumed<E: SolarGlobalEphemeris + ?Sized>(
    eph: &E,
    jd_new_moon: f64,
    mask: i32,
) -> Result<Option<(i32, [f64; 10])>> {
    let moon_lon = eph.candidate_moon_longitude_deg(jd_new_moon)?;
    let node_lon = eph.mean_node_lon_deg(jd_new_moon)?;
    let node_dist = moon_node_distance_deg(node_lon, moon_lon)?;
    if node_dist < SOLAR_NODE_PREFILTER_DEG {
        let jd_max = {
            let mut gamma_at = |jd: f64| eph.gamma_er(jd);
            refine_solar_eclipse_maximum(
                &mut gamma_at,
                jd_new_moon,
                SOLAR_MAXIMUM_SEARCH_RANGE_DAYS,
            )?
        };
        let state_max = eph.sun_moon_state(jd_max)?;
        let gamma_max = eph.gamma_er(jd_max)?;
        let (l1, l2) = eph.shadow_limits_er(jd_max)?;
        let mut where_at_max = |magnitude| {
            if magnitude {
                eph.where_magnitude(jd_max).map(|value| (0, value))
            } else {
                eph.where_retflag(jd_max).map(|value| (value, 0.))
            }
        };
        let mut hybrid_at_max = || {
            let mut contact = |jd: f64, target: f64, before: bool| -> Result<f64> {
                let mut gamma_at = |jd: f64| eph.gamma_er(jd);
                find_contact_time_gamma(
                    &mut gamma_at,
                    jd,
                    target,
                    before,
                    HYBRID_CONTACT_SEARCH_RANGE_DAYS,
                )
            };
            let mut width = |jd: f64| eph.ground_core_width_km(jd);
            is_hybrid_solar_eclipse(jd_max, l2, &mut contact, &mut width)
        };
        let classified = classify_solar_eclipse_consumed(
            &state_max,
            gamma_max,
            l1,
            l2,
            &mut where_at_max,
            &mut hybrid_at_max,
        )?;
        if classified.eclipse_type != 0 && sol_glob_accepts(mask, classified.eclipse_type) {
            let mut residuals = |jd: f64| eph.contact_residuals_km(jd);
            let mut x = |jd: f64| eph.besselian_x_er(jd);
            let times = calculate_solar_phases_besselian(
                jd_max,
                classified.eclipse_type,
                &mut residuals,
                &mut x,
            )?;
            return Ok(Some((classified.eclipse_type, times)));
        }
    }
    Ok(None)
}

#[cfg(test)]
mod refined_disc_contact_tests {
    use super::*;

    #[test]
    fn shallow_entry_and_exit_contacts_choose_the_best_representable_epoch() {
        let origin: f64 = 2460000.0;
        let target = 0.000003141592653589793;
        for increasing in [false, true] {
            let sign = if increasing { 1.0 } else { -1.0 };
            let slope = sign * 1e-4;
            let exact = origin + sign * target / slope;
            let mut separation = |t| Ok((t - origin) * slope);
            let actual = details_contact_time(
                &mut separation,
                origin,
                origin + 0.1,
                sign * target,
                increasing,
            )
            .unwrap();
            let residual = |t: f64| ((t - origin) * slope - sign * target).abs();
            assert!((actual - exact).abs() <= exact.next_up() - exact);
            assert!(residual(actual) <= residual(actual.next_down()));
            assert!(residual(actual) <= residual(actual.next_up()));
        }
    }

    #[test]
    fn already_adjacent_contacts_select_an_actual_endpoint() {
        let lo: f64 = 2460000.0;
        let hi = lo.next_up();
        let mut separation = |t| Ok(t - lo);
        let actual = details_contact_time(&mut separation, lo, hi, 0.75 * (hi - lo), true).unwrap();
        assert_eq!(actual, hi);
    }

    #[test]
    fn nonfinite_contact_geometry_is_an_error() {
        assert!(matches!(
            details_contact_time(&mut |_| Ok(f64::NAN), 0.0, 1.0, 0.5, true),
            Err(Error::NonFiniteSearchValue { .. })
        ));
    }
}
