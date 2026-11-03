//! Output-frame selection, reported rates and coordinate shaping for stars.
//!
//! This layer consumes an already acquired six-slot place. Its borrowed serving
//! adapter supplies the existing frame/model/clock inputs only when requested by
//! a flag branch; row, observer and raw-place acquisition stay with serving.

use super::motion::sin_cos;
use super::{
    Error, FixstarFlagEnvironment, STAR_FLAG_LAYER_STEP_DAYS, STAR_J2000_OBLIQUITY_DEG,
    SealedStarServing,
};

/// Backward precession of ecliptic coordinates between two epochs
///
/// 1. ecliptic→equatorial at the from-epoch mean obliquity;
/// 2. equatorial precession `P(to) @ P(from)^T` over the native
///    bias-free Vondrák matrices;
/// 3. equatorial→ecliptic at the to-epoch mean obliquity with
///    `atan2`/`asin` shaping and `[0, 360)` normalization.
///
/// Returns the input unchanged when `|jd_tt - J2000| < 1e-6`.
/// Errors: `Error::NonFiniteArgument` on non-finite input; frame and shaping errors propagate from the obliquity, precession-matrix, and rotation calls.
pub fn precess_ecliptic_to_j2000(
    lon_deg: f64,
    lat_deg: f64,
    jd_tt: f64,
) -> Result<(f64, f64), Error> {
    use crate::constants::J2000;
    use crate::positions::{
        mat3_mul, mat3_vec3, rotate_equatorial_to_ecliptic, vondrak_mean_obliquity_deg,
        vondrak_precession_matrix_nobias,
    };
    for v in [lon_deg, lat_deg, jd_tt] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    if (jd_tt - J2000).abs() < 1e-6 {
        return Ok((lon_deg, lat_deg));
    }
    let eps1 = vondrak_mean_obliquity_deg(jd_tt)?.to_radians();
    let eps2 = vondrak_mean_obliquity_deg(J2000)?.to_radians();
    // Step 1: ecliptic→equatorial at the from-epoch mean obliquity
    let lon_rad = lon_deg.to_radians();
    let lat_rad = lat_deg.to_radians();
    let (sin_lon, cos_lon) = sin_cos(lon_rad);
    let (sin_lat, cos_lat) = sin_cos(lat_rad);
    let (sin_e1, cos_e1) = sin_cos(eps1);
    let x1 = cos_lat * cos_lon;
    let y1 = cos_lat * sin_lon * cos_e1 - sin_lat * sin_e1;
    let z1 = cos_lat * sin_lon * sin_e1 + sin_lat * cos_e1;
    // Step 2: Vondrák 2011 equatorial precession from `jd_tt` to J2000.
    let p_to = vondrak_precession_matrix_nobias(J2000)?;
    let p_from = vondrak_precession_matrix_nobias(jd_tt)?;
    let p_from_t = [
        [p_from[0][0], p_from[1][0], p_from[2][0]],
        [p_from[0][1], p_from[1][1], p_from[2][1]],
        [p_from[0][2], p_from[1][2], p_from[2][2]],
    ];
    let prec = mat3_mul(p_to, p_from_t)?;
    let eq2 = mat3_vec3(prec, [x1, y1, z1])?;
    // Step 3: equatorial→ecliptic at the J2000 mean obliquity.
    let ecl = rotate_equatorial_to_ecliptic(eq2, eps2)?;
    let dist = (ecl[0] * ecl[0] + ecl[1] * ecl[1] + ecl[2] * ecl[2]).sqrt();
    let new_lon = ecl[1].atan2(ecl[0]).to_degrees().rem_euclid(360.0);
    let new_lat = if dist > 0.0 {
        (ecl[2] / dist).clamp(-1.0, 1.0).asin().to_degrees()
    } else {
        0.0
    };
    Ok((new_lon, new_lat))
}

/// Reduce a single-epoch ecliptic star position through coordinate/frame flags.
///
/// 1. J2000 backward precession, else NONUT/SIDEREAL mean-ecliptic
///    reduction via nutation-in-longitude subtraction;
/// 2. EQUATORIAL rotation with the J2000 / mean / true obliquity arm;
/// 3. SIDEREAL mean-ayanamsha subtraction from the first coordinate
///    (longitude, or right ascension after the equatorial rotation).
///
/// Longitude/latitude and obliquities are degrees; distance is AU and passes
/// through unchanged. The equatorial rotation uses the direct expression
/// `sin(lat') = sin(lat)*cos(eps) + cos(lat)*sin(eps)*sin(lon)` and full
/// Cartesian longitude components, without normalizing an intermediate
/// vector. J2000 takes precedence over NONUT/SIDEREAL; its obliquity is
/// 84381.406 arcseconds. Epoch `jd_tt` is TT Julian Day. Non-finite inputs or
/// obliquity return [`Error::NonFiniteArgument`]; environment errors propagate.
pub fn fixstar_positional_transform_via(
    lon_deg: f64,
    lat_deg: f64,
    dist_au: f64,
    jd_tt: f64,
    flags: i32,
    j2000_native: bool,
    env: &impl FixstarFlagEnvironment,
) -> Result<[f64; 3], Error> {
    use crate::constants::{FLG_EQUATORIAL, FLG_J2000, FLG_NONUT, FLG_SIDEREAL};
    for v in [lon_deg, lat_deg, dist_au, jd_tt] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    let (mut plon, mut plat) = (lon_deg, lat_deg);
    // Step 1: frame selection (J2000 takes precedence over the
    // NONUT/SIDEREAL branch, `if/elif`).
    if flags & FLG_J2000 != 0 {
        if !j2000_native {
            (plon, plat) = env.precess_ecliptic_to_j2000(plon, plat, jd_tt)?;
        }
    } else if flags & (FLG_NONUT | FLG_SIDEREAL) != 0 {
        let dpsi_rad = env.nutation_dpsi_rad(jd_tt)?;
        plon = (plon - dpsi_rad.to_degrees()).rem_euclid(360.0);
    }
    // Step 2: equatorial rotation.
    if flags & FLG_EQUATORIAL != 0 {
        let eps = if flags & FLG_J2000 != 0 {
            STAR_J2000_OBLIQUITY_DEG
        } else if flags & (FLG_NONUT | FLG_SIDEREAL) != 0 {
            env.mean_obliquity_deg(jd_tt)?
        } else {
            env.true_obliquity_deg(jd_tt)?
        };
        if !eps.is_finite() {
            return Err(Error::NonFiniteArgument { value: eps });
        }
        // The public spherical rotation retains its direct sine-latitude
        // expression; the fast-calc normalized-Cartesian route is distinct.
        let (lon, lat, _) = crate::events::cotrans(plon, plat, dist_au, -eps);
        // Canonicalize a wrapped zero angle to positive zero.
        (plon, plat) = (if lon == 0.0 { 0.0 } else { lon }, lat);
    }
    // Step 3: sidereal mean-ayanamsha subtraction (deliberately NOT
    // renormalized past a zero crossing: the result comes from `rem_euclid`,
    // and the mean (not true) value is removed).
    if flags & FLG_SIDEREAL != 0 {
        use crate::constants::{FLG_NOABERR, FLG_NOGDEFL, FLG_TRUEPOS};
        let jd_ut = env.ut1_for_tt(jd_tt)?;
        let ab_bits = flags & (FLG_TRUEPOS | FLG_NOABERR | FLG_NOGDEFL);
        let mean_ayan = env.mean_ayanamsha_deg_ut(jd_ut, FLG_NONUT | ab_bits)?;
        plon = (plon - mean_ayan).rem_euclid(360.0);
    }
    Ok([plon, plat, dist_au])
}

/// Full post-calculation flag application: reported position at the request epoch, angular speeds by
/// central-differencing the reported angles (raw endpoints reconstructed
/// as midpoint ± raw_speed · `h` with the fixed `h = 0.5`, each reduced at
/// its own epoch so obliquity/nutation/ayanamsha drift rates survive —
/// `rem_euclid` shortest-arc fold on the first coordinate), then XYZ
/// (priority) or RADIANS output shaping. The geometric radial derivative is
/// invariant under these angle transforms and passes through unchanged.
/// Errors: `Error::NonFiniteArgument` on non-finite raw place or epoch; per-epoch transform and Cartesian-shaping errors propagate.
pub fn apply_fixstar_flags_via(
    raw: [f64; 6],
    jd_tt: f64,
    flags: i32,
    j2000_native: bool,
    env: &impl FixstarFlagEnvironment,
) -> Result<[f64; 6], Error> {
    use crate::constants::{FLG_RADIANS, FLG_SPEED, FLG_XYZ};
    for v in raw {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let center =
        fixstar_positional_transform_via(raw[0], raw[1], raw[2], jd_tt, flags, j2000_native, env)?;
    let (speed_lon, speed_lat, speed_dist) = if flags & FLG_SPEED != 0 {
        let h = STAR_FLAG_LAYER_STEP_DAYS;
        let plus = fixstar_positional_transform_via(
            raw[0] + raw[3] * h,
            raw[1] + raw[4] * h,
            raw[2] + raw[5] * h,
            jd_tt + h,
            flags,
            j2000_native,
            env,
        )?;
        let minus = fixstar_positional_transform_via(
            raw[0] - raw[3] * h,
            raw[1] - raw[4] * h,
            raw[2] - raw[5] * h,
            jd_tt - h,
            flags,
            j2000_native,
            env,
        )?;
        let d_lon = (plus[0] - minus[0] + 180.0).rem_euclid(360.0) - 180.0;
        (
            d_lon / (2.0 * h),
            (plus[1] - minus[1]) / (2.0 * h),
            // These frame/sidereal transforms change angles only. Preserve
            // the physical rate instead of rounding it into million-AU
            // endpoints and subtracting those rounded distances again.
            raw[5],
        )
    } else {
        (0.0, 0.0, 0.0)
    };
    if flags & FLG_XYZ != 0 {
        return Ok(crate::positions::spherical_to_cartesian_with_velocity(
            center[0], center[1], center[2], speed_lon, speed_lat, speed_dist,
        ));
    }
    if flags & FLG_RADIANS != 0 {
        return Ok([
            center[0].to_radians(),
            center[1].to_radians(),
            center[2],
            speed_lon.to_radians(),
            speed_lat.to_radians(),
            speed_dist,
        ]);
    }
    Ok([
        center[0], center[1], center[2], speed_lon, speed_lat, speed_dist,
    ])
}

/// Return flags echoed by the star entries: UT entries resolve center priority, force exactly one
/// ephemeris bit and add the implied correction bits; TT entries echo the
/// request.
pub fn fixstar_ret_flags(flags_in: i32, implied: bool) -> i32 {
    use crate::positions::{exclusive_ephemeris_bit, implied_retflag_bits, resolve_center_flags};
    let mut ret = flags_in;
    if implied {
        ret = resolve_center_flags(ret);
        ret = exclusive_ephemeris_bit(ret, false);
        ret |= implied_retflag_bits(ret);
    }
    ret
}

/// Borrowed frame/model environment for the public star flag layer.
///
/// Mean obliquity is Vondrak 2011; true obliquity adds native IAU 2006/2000 A
/// nutation, independently of stored LEB nutation. TT-to-UT1 uses the serving's
/// explicit bundled or callback clock. Mean ayanamsha reacquires evaluation TT
/// through that same star clock, with the bound dynamic-anchor provider.
/// Model Delta-T supplies defining epochs separately. No provider borrow
/// escapes, and consumed frame, clock, anchor or model errors propagate.
pub(super) struct SealedFlagEnv<'b, P: crate::positions::StateProvider> {
    // Parent serving composes this view only after acquiring the raw place.
    pub(super) serving: &'b SealedStarServing<'b, P>,
}

impl<P: crate::positions::StateProvider> FixstarFlagEnvironment for SealedFlagEnv<'_, P> {
    fn precess_ecliptic_to_j2000(
        &self,
        lon_deg: f64,
        lat_deg: f64,
        jd_tt: f64,
    ) -> Result<(f64, f64), Error> {
        precess_ecliptic_to_j2000(lon_deg, lat_deg, jd_tt)
    }

    fn nutation_dpsi_rad(&self, jd_tt: f64) -> Result<f64, Error> {
        let (dpsi, _) = crate::erfa::nut06a(jd_tt)?;
        Ok(dpsi)
    }

    fn mean_obliquity_deg(&self, jd_tt: f64) -> Result<f64, Error> {
        crate::positions::vondrak_mean_obliquity_deg(jd_tt)
    }

    fn true_obliquity_deg(&self, jd_tt: f64) -> Result<f64, Error> {
        let (_, deps) = crate::erfa::nut06a(jd_tt)?;
        let mean = crate::positions::vondrak_mean_obliquity_deg(jd_tt)?;
        Ok(mean + deps.to_degrees())
    }

    fn ut1_for_tt(&self, jd_tt: f64) -> Result<f64, Error> {
        self.serving.ut1_for_tt(jd_tt)
    }

    fn mean_ayanamsha_deg_ut(&self, jd_ut: f64, ex_flags: i32) -> Result<f64, Error> {
        use crate::constants::{FLG_NOABERR, FLG_NOGDEFL, FLG_TRUEPOS};
        let anchors = self.serving.anchor_provider();
        let aya = crate::sidereal::SealedAyanamshaServing::new(
            self.serving.sidereal,
            self.serving.deltat_days_at,
            Some(&anchors),
        );
        // Evaluation TT through the serving clock: routing this arm through
        // the model Delta-T callback instead shifts Fagan-default
        // longitudes by ~2.1 e-9 deg at modern epochs. Anchor toggles
        // follow the `ex_ut` contract (flags alone, no helio/bary
        // centering); `FLG_NONUT` always selects the mean value here.
        let jd_tt = self.serving.jd_tt_for_ut(jd_ut)?;
        let noaberr = ex_flags & (FLG_TRUEPOS | FLG_NOABERR) != 0;
        let nogdefl = ex_flags & (FLG_TRUEPOS | FLG_NOGDEFL) != 0;
        let (mean, _) = aya.mean_tt_with(jd_tt, noaberr, nogdefl)?;
        Ok(mean)
    }
}
