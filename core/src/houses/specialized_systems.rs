//! Specialized systems implementation for the parent module.
//! Public entry points and types are exposed through the parent API.

use super::*;

/// Natural-gradient houses: Equal with 0° Aries as cusp 1 — the Whole-Sign-from-Aries
/// reading — with `λᵢ = (i−1)·30°` (H1 = 0° through H12 = 330°).
/// Time/location independent; works at all latitudes; exact; all five
/// inputs only shape the call. Returns `[f64; 12]` houses 1–12 (the
/// release 13-list with index 0 unused, rendered in Rust order).
pub fn houses_natural_gradient(
    _armc_deg: f64,
    _lat_deg: f64,
    _eps_deg: f64,
    _asc_deg: f64,
    _mc_deg: f64,
) -> [f64; 12] {
    let mut cusps = [0.0_f64; 13];
    for (n, cusp) in cusps.iter_mut().enumerate().skip(1) {
        *cusp = wrap360((n as f64 - 1.0) * 30.0);
    }
    strip_index0(cusps)
}

/// Horizontal (azimuthal) houses: the horizon divided into twelve 30° arcs from the east point,
/// each division carried to the ecliptic along its vertical circle
/// (Holden, *Elements of House Division*, ch. 10; rising-point
/// derivation after Smart ch. 3).
///
/// Signed co-latitude `φ' = ±90° − φ` with the equatorial `|φ'| = 90°`
/// singularity clamped; each division at `θ = 30°/60°/90°/120°/150°`
/// maps to cusps 3/2/1/12/11 with pole height
/// `arcsin(sin θ·sin φ')` and RA offset
/// `atan2(cos θ, sin θ·cos φ')` from the `ARMC + 270°` base (`ARMC +
/// 180°` with the `+90°` call convention), evaluated through
/// [`calc_ascendant`]. The polar-circle branch (Asc/MC flip), the 180°
/// flip of the computed cusps 1/2/3/11/12, the Asc/DC orientation check
/// and the MC/IC + opposite derivation follow the release.
///
/// Arguments are degrees (active ARMC, latitude, true obliquity,
/// precomputed Ascendant/MC); takes the passed angles as-is (entry owns
/// the R/C/T/J MC-flip and ASCMC assembly). Returns `[f64; 12]` houses
/// 1–12.
pub fn houses_horizontal(
    armc_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
    asc_deg: f64,
    mc_deg: f64,
) -> [f64; 12] {
    let mut cusps = [0.0_f64; 13];
    let mut asc = asc_deg;
    let mut mc = mc_deg;

    let mut co_lat = if lat_deg > 0.0 {
        90.0 - lat_deg
    } else {
        -90.0 - lat_deg
    };
    if (co_lat.abs() - 90.0).abs() < HORIZONTAL_NEAR_ZERO_DEG {
        co_lat = if co_lat > 0.0 {
            90.0 - HORIZONTAL_NEAR_ZERO_DEG
        } else {
            -90.0 + HORIZONTAL_NEAR_ZERO_DEG
        };
    }

    let armc_base = wrap360(armc_deg + 180.0);
    let co_lat_rad = co_lat.to_radians();
    let (sin_colat, cos_colat) = sin_cos(co_lat_rad);

    const CUSP_MAP: [(usize, f64); 5] = [(3, 30.0), (2, 60.0), (1, 90.0), (12, 120.0), (11, 150.0)];
    for (cusp_idx, theta_deg) in CUSP_MAP {
        let theta_rad = theta_deg.to_radians();
        let (sin_theta, cos_theta) = sin_cos(theta_rad);
        let pole_height = clamp_unit(sin_theta * sin_colat).asin().to_degrees();
        let ra_offset = cos_theta.atan2(sin_theta * cos_colat).to_degrees();
        cusps[cusp_idx] =
            calc_ascendant(armc_base + 90.0 + ra_offset, eps_deg, lat_deg, pole_height);
    }

    if co_lat.abs() >= 90.0 - eps_deg {
        let asc_mc_diff = wrap360(asc - mc + 540.0) - 180.0;
        if asc_mc_diff < 0.0 {
            asc = wrap360(asc + 180.0);
            mc = wrap360(mc + 180.0);
            for (n, cusp) in cusps.iter_mut().enumerate().skip(1) {
                if (4..10).contains(&n) {
                    continue;
                }
                *cusp = wrap360(*cusp + 180.0);
            }
        }
    }

    for i in [1_usize, 2, 3, 11, 12] {
        cusps[i] = wrap360(cusps[i] + 180.0);
    }

    let asc_mc_diff = wrap360(asc - mc + 540.0) - 180.0;
    if asc_mc_diff < 0.0 {
        asc = wrap360(asc + 180.0);
    }
    // No later kernel statement reads `asc` (house 7 derives from cusp 1
    // below); the trailing read keeps the release reassignment above
    // without a dead store. The release has the same shape.
    let _ = asc;

    cusps[10] = mc;
    cusps[4] = wrap360(mc + 180.0);
    cusps[7] = wrap360(cusps[1] + 180.0);
    cusps[8] = wrap360(cusps[2] + 180.0);
    cusps[9] = wrap360(cusps[3] + 180.0);
    cusps[5] = wrap360(cusps[11] + 180.0);
    cusps[6] = wrap360(cusps[12] + 180.0);
    strip_index0(cusps)
}

/// One APC house cusp.
///
/// The APC construction is the Sunshine construction applied to the
/// Ascendant's own parallel of declination: the parallel's diurnal and
/// nocturnal arcs are divided into six equal parts each and the division
/// points are carried to the ecliptic along house circles through the
/// north/south horizon axis (de Boer, "APC Houses", ¶1–3). The eastern
/// half selects the Ascendant vector (the raw closed form can pick the
/// setting point inside the polar circles); the hour angle wraps to
/// `(−π, π]`; the day-semiarc rule and the trisection table give the
/// anchor, with the four explicit angle assignments. Inputs are radians
/// (the caller works in radians); output is degrees in `[0, 360)`.
/// Cusp 1 reproduces the Ascendant, and cusp 10's house circle is the
/// meridian by construction: at index 10 the hour angle is 0, so the
/// anchor transits the upper meridian and the north/south circle
/// through it is the meridian itself. Scalar on both sides.
pub fn apc_cusp(index: usize, lat_rad: f64, eps_rad: f64, armc_rad: f64) -> f64 {
    let (zenith, east, north) = horizon_frame(armc_rad, lat_rad);
    let ecl_pole: EqVec = {
        let (sin_eps, cos_eps) = sin_cos(eps_rad);
        [0.0, -sin_eps, cos_eps]
    };
    let mut asc_vec = vec_normalized(vec_cross(ecl_pole, zenith));
    if vec_dot(asc_vec, east) < 0.0 {
        asc_vec = [-asc_vec[0], -asc_vec[1], -asc_vec[2]];
    }
    let dec_asc = clamp_unit(asc_vec[2]).asin();
    let ra_asc = asc_vec[1].atan2(asc_vec[0]);
    let ha_asc = {
        let diff = armc_rad - ra_asc;
        let (sin_diff, cos_diff) = sin_cos(diff);
        sin_diff.atan2(cos_diff)
    };
    let day_semiarc = if ha_asc < 0.0 {
        -ha_asc
    } else {
        2.0 * std::f64::consts::PI - ha_asc
    };

    let mut hour_angle = 0.0_f64;
    for (idx, angle) in semiarc_cusp_hour_angles(day_semiarc) {
        if idx == index {
            hour_angle = angle;
            break;
        }
    }
    // The four angles bypass the trisection table ( assignments).
    match index {
        1 => hour_angle = -day_semiarc,
        7 => hour_angle = day_semiarc,
        10 => hour_angle = 0.0,
        4 => hour_angle = std::f64::consts::PI,
        _ => {}
    }
    let anchor = unit_radec(armc_rad - hour_angle, dec_asc);
    house_circle_lon_deg(anchor, north, eps_rad, true)
}

/// APC cusp with a checked one-based cusp index.
///
/// The index must be in `1..=12`; latitude, obliquity and ARMC are radians.
/// Returns ecliptic longitude in degrees using [`apc_cusp`]'s native geometry.
/// An invalid index returns [`Error::InvalidHouseCuspIndex`]. Degenerate
/// geometric inputs retain the underlying construction's binary64 behavior.
pub fn apc_cusp_checked(index: i64, lat_rad: f64, eps_rad: f64, armc_rad: f64) -> Result<f64> {
    if !(1..=12).contains(&index) {
        return Err(Error::InvalidHouseCuspIndex { index });
    }
    Ok(apc_cusp(index as usize, lat_rad, eps_rad, armc_rad))
}

/// APC houses: the Ascendant-parallel
/// 6+6 division with north/south position circles (Knegt,
/// *Handleiding*): every cusp from [`apc_cusp`], then the directly
/// computed MC substituted for cusp 10 (the per-cusp formula is
/// ill-conditioned near `|lat| = 90°`) with cusp 4 its opposite, and the
/// polar-circle horizon-crossing flip (`difdeg2n(asc, mc) < 0` → all
/// cusps + 180°).
///
/// Arguments are degrees (active ARMC, latitude, true obliquity,
/// precomputed Ascendant/MC); takes the passed angles as-is (the APC
/// `ascmc[1]` re-sync stays entry scope). Returns `[f64; 12]` houses
/// 1–12.
pub fn houses_apc(
    armc_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
    asc_deg: f64,
    mc_deg: f64,
) -> [f64; 12] {
    let lat_rad = lat_deg.to_radians();
    let eps_rad = eps_deg.to_radians();
    let armc_rad = armc_deg.to_radians();
    let mut cusps = [0.0_f64; 13];
    for (n, cusp) in cusps.iter_mut().enumerate().skip(1) {
        *cusp = apc_cusp(n, lat_rad, eps_rad, armc_rad);
    }
    cusps[10] = mc_deg;
    cusps[4] = wrap360(mc_deg + 180.0);

    if lat_deg.abs() >= 90.0 - eps_deg && difdeg2n_deg(asc_deg, mc_deg) < 0.0 {
        // The release also reassigns its locals `asc`/`mc` here; only the
        // flipped MC reaches the caller, through `ascmc` ( entry
        // scope). The kernel output is the flipped cusp array below.
        for cusp in cusps.iter_mut().skip(1) {
            *cusp = wrap360(*cusp + 180.0);
        }
    }
    strip_index0(cusps)
}

/// Krusinski-Pisa houses: the Ascendant–zenith great circle in twelve 30° arcs from
/// the Ascendant toward the zenith (Krusinski's own 1995 definition),
/// the zenith arc at 90° giving cusp 10, each division carried to the
/// ecliptic along its hour circle. `mc` only shapes the call — cusp 10
/// is the zenith's division point, constructed, never substituted.
/// Arguments are degrees (active ARMC, latitude, true obliquity,
/// Ascendant); takes the passed angles as-is. Returns `[f64; 12]`
/// houses 1–12.
/// Returns [`Error::ZeroVectorNormalize`] when the Ascendant direction
/// coincides with the zenith (the Ascendant–zenith circle is
/// undefined): the release raises `ZeroDivisionError` from the same
/// normalization step, never quiet NaNs.
pub fn houses_krusinski(
    armc_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
    asc_deg: f64,
    _mc_deg: f64,
) -> Result<[f64; 12]> {
    let eps_rad = eps_deg.to_radians();
    let asc_rad = asc_deg.to_radians();
    let (sin_asc, cos_asc) = sin_cos(asc_rad);
    let (sin_eps, cos_eps) = sin_cos(eps_rad);
    let a_vec: EqVec = [cos_asc, sin_asc * cos_eps, sin_asc * sin_eps];
    let zenith = unit_radec(armc_deg.to_radians(), lat_deg.to_radians());
    // Component of the zenith orthogonal to the Ascendant: second basis
    // vector of the Ascendant-zenith circle. The Ascendant lies on the
    // horizon, so the two are already orthogonal up to rounding.
    let z_dot_a = vec_dot(zenith, a_vec);
    let orth: EqVec = [
        zenith[0] - z_dot_a * a_vec[0],
        zenith[1] - z_dot_a * a_vec[1],
        zenith[2] - z_dot_a * a_vec[2],
    ];
    //  release gate: `_normalized` divides by
    // `sqrt(v₀² + v₁² + v₂²)`; a zero norm has no defined direction.
    // exactly when that norm is `0.0` (any numerator). The squared-norm
    // spelling below is the same gate (`n == 0 ⟺ n² == 0`).
    if orth[0] * orth[0] + orth[1] * orth[1] + orth[2] * orth[2] == 0.0 {
        return Err(Error::ZeroVectorNormalize);
    }
    let up = vec_normalized(orth);
    let mut cusps = [0.0_f64; 13];
    const ARCS: [(usize, f64); 6] = [
        (1, 0.0),
        (12, 30.0),
        (11, 60.0),
        (10, 90.0),
        (2, -30.0),
        (3, -60.0),
    ];
    for (cusp, arc_deg) in ARCS {
        let arc = arc_deg.to_radians();
        let (sin_arc, cos_arc) = sin_cos(arc);
        let point: EqVec = [
            cos_arc * a_vec[0] + sin_arc * up[0],
            cos_arc * a_vec[1] + sin_arc * up[1],
            cos_arc * a_vec[2] + sin_arc * up[2],
        ];
        cusps[cusp] = hour_circle_lon_deg(point, eps_rad);
    }
    for (src, opp) in [
        (1_usize, 7_usize),
        (2, 8),
        (3, 9),
        (10, 4),
        (11, 5),
        (12, 6),
    ] {
        cusps[opp] = wrap360(cusps[src] + 180.0);
    }
    Ok(strip_index0(cusps))
}

/// Savard-A houses: the prime vertical from the east point (declination 0) to the zenith
/// (declination = latitude); the declination parallels at one and two
/// thirds of the latitude anchor position circles through the horizon
/// north/south points (Savard, "Astrological House Systems",
/// "Albategnus" entry). Two-thirds anchor → cusps 11/5, one-third →
/// 12/6, west-side mirrors → 3/9 and 2/8; the anchor at declination
/// `k·lat/3` sits at prime-vertical arc
/// `sin(t) = sin(k·lat/3)/sin(lat)` (equatorial limit `k/3`). At polar
/// latitudes the east/west numbering ties to the actual rising point,
/// not to the horizon frame's representational choice. Arguments are
/// degrees (active ARMC, latitude, true obliquity, Ascendant = cusp 1,
/// MC = cusp 10); takes the passed angles as-is. Returns `[f64; 12]`
/// houses 1–12.
pub fn houses_savard_a(
    armc_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
    asc_deg: f64,
    mc_deg: f64,
) -> [f64; 12] {
    let eps_rad = eps_deg.to_radians();
    let lat_rad = lat_deg.to_radians();
    let (zenith, east, north) = horizon_frame(armc_deg.to_radians(), lat_rad);
    let asc_rad = asc_deg.to_radians();
    let (sin_asc, cos_asc) = sin_cos(asc_rad);
    let (sin_eps, cos_eps) = sin_cos(eps_rad);
    let asc_vec: EqVec = [cos_asc, sin_asc * cos_eps, sin_asc * sin_eps];
    let frame_east_sign = if vec_dot(east, asc_vec) >= 0.0 {
        1.0
    } else {
        -1.0
    };

    let mut cusps = [0.0_f64; 13];
    cusps[1] = asc_deg;
    cusps[10] = mc_deg;
    let sin_lat = lat_rad.sin();
    for (third, cusp_east, cusp_west) in [(2_i32, 11_usize, 3_usize), (1, 12, 2)] {
        let s = if sin_lat.abs() < 1e-10 {
            third as f64 / 3.0
        } else {
            (third as f64 * lat_deg / 3.0).to_radians().sin() / sin_lat
        };
        let t = clamp_unit(s).asin();
        let (sin_t, cos_t) = sin_cos(t);
        for (east_sign, cusp) in [(1.0_f64, cusp_east), (-1.0, cusp_west)] {
            let anchor: EqVec = [
                east_sign * cos_t * east[0] + sin_t * zenith[0],
                east_sign * cos_t * east[1] + sin_t * zenith[1],
                east_sign * cos_t * east[2] + sin_t * zenith[2],
            ];
            cusps[cusp] =
                house_circle_lon_deg(anchor, north, eps_rad, east_sign * frame_east_sign > 0.0);
        }
    }
    for (src, opp) in [
        (1_usize, 7_usize),
        (10, 4),
        (11, 5),
        (12, 6),
        (2, 8),
        (3, 9),
    ] {
        cusps[opp] = wrap360(cusps[src] + 180.0);
    }
    strip_index0(cusps)
}

/// One trisected semi-arc offset mapped to an ecliptic longitude
/// the 5-step spherical-triangle
/// method (Makransky 1988, *Primary Directions*): great-circle arc
/// `xhs = 2·arcsin(cos δ·sin(offset/2))`;
/// `cos α = tan δ·tan(xhs/2)`; the diurnal vs nocturnal triangle
/// side/angle; cosine-law side `c`; sine-law zenith distance with the
/// `1e-6` floor; RA intersection with signed pole height;
/// rising-point projection via [`ra_to_ecliptic_longitude`]. Diurnal
/// offsets map to the houses 8–12 side, nocturnal offsets to the
/// houses 2–6 side.
///
/// The release passes the pre-computed sines as keyword arguments; the
/// Rust form takes them positionally in the same order (arc offset and
/// ARMC/latitude/Sun-declination in degrees, flags and trig legs as
/// given). Returns degrees in `[0, 360)`; scalar on both sides.
#[allow(
    clippy::too_many_arguments,
    reason = "verbatim release `_sunshine_arc_to_ecliptic` call shape \
    (arc offset, diurnal flag, armc/lat/sun_dec plus the six pre-computed \
    sines); grouping would hide the operation-order audit"
)]
pub(super) fn sunshine_arc_to_ecliptic(
    arc_offset_ra_deg: f64,
    is_diurnal: bool,
    armc_deg: f64,
    lat_deg: f64,
    sun_dec_deg: f64,
    sin_lat: f64,
    cos_lat: f64,
    cos_sun_dec: f64,
    tan_sun_dec: f64,
    sin_obliquity: f64,
    cos_obliquity: f64,
) -> f64 {
    let great_circle_arc = 2.0
        * (cos_sun_dec * (arc_offset_ra_deg / 2.0).to_radians().sin())
            .asin()
            .to_degrees();
    let cos_alpha = clamp_unit(tan_sun_dec * (great_circle_arc / 2.0).to_radians().tan());
    let alpha = cos_alpha.acos().to_degrees();
    let (triangle_angle, triangle_side) = if is_diurnal {
        (180.0 - alpha, 90.0 - lat_deg + sun_dec_deg)
    } else {
        (alpha, 90.0 - lat_deg - sun_dec_deg)
    };
    let great_rad = great_circle_arc.to_radians();
    let (sin_great, cos_great) = sin_cos(great_rad);
    let side_rad = triangle_side.to_radians();
    let (sin_side, cos_side) = sin_cos(side_rad);
    let angle_rad = triangle_angle.to_radians();
    let (sin_angle, cos_angle) = sin_cos(angle_rad);
    let cos_c = clamp_unit(cos_great * cos_side + sin_great * sin_side * cos_angle);
    let side_c = cos_c.acos().to_degrees();
    let zenith_dist = if side_c < 1e-6 {
        0.0
    } else {
        clamp_unit(sin_great * sin_angle / side_c.to_radians().sin())
            .asin()
            .to_degrees()
    };
    let equator_ra_offset = (cos_lat * zenith_dist.to_radians().tan())
        .atan()
        .to_degrees();
    let mut pole_height = (zenith_dist.to_radians().sin() * sin_lat)
        .asin()
        .to_degrees();
    let ra_intersection = if is_diurnal {
        wrap360(armc_deg + equator_ra_offset)
    } else {
        pole_height = -pole_height;
        wrap360(equator_ra_offset + armc_deg + 180.0)
    };
    ra_to_ecliptic_longitude(ra_intersection, pole_height, sin_obliquity, cos_obliquity)
}

/// Sunshine/Treindl houses: the
/// Makransky public Sunshine definition realized by the Treindl-style
/// numerical construction (code `I` here; code `i` selects the distinct
/// Makransky upper/lower-meridian construction in [`houses_sunshine_makransky`]):
/// the Sun's
/// diurnal semi-arc (`DSA = 90° + AD`) and nocturnal semi-arc
/// (`NSA = 90° − AD`) from the ascensional difference
/// (`sin AD = tan δ☉·tan φ`, clamped at `±(90° − 1e-10)`), trisected to
/// RA offsets and mapped through `sunshine_arc_to_ecliptic`
/// (nocturnal offsets → houses 2/3/5/6, diurnal → 8/9/11/12). Cusps
/// 11/12/2/3 are generally NOT opposite 5/6/8/9. The MC-below-horizon
/// flip (180° on MC/IC and the eight intermediate cusps, decided by the
/// `|lat − mcdec| > 90°` gate) selects the hemisphere of the intermediate arcs.
///
/// `sun_dec_deg` is the explicit Sun declination; remaining arguments are degrees
/// (active ARMC, latitude, true obliquity, precomputed Ascendant/MC).
/// Takes the passed angles as-is. Returns `[f64; 12]` houses 1–12 (the
/// release 13-list with index 0 unused, rendered in Rust order).
pub fn houses_sunshine(
    armc_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
    asc_deg: f64,
    mc_deg: f64,
    sun_dec_deg: f64,
) -> [f64; 12] {
    let mcdec = (armc_deg.to_radians().sin() * eps_deg.to_radians().tan())
        .atan()
        .to_degrees();
    let mc_under_horizon = (lat_deg - mcdec).abs() > 90.0;

    let mut mc = mc_deg;
    if mc_under_horizon {
        mc = wrap360(mc + 180.0);
    }
    let mut cusps = [0.0_f64; 13];
    cusps[1] = asc_deg;
    cusps[10] = mc;
    cusps[4] = wrap360(mc + 180.0);
    cusps[7] = wrap360(asc_deg + 180.0);

    let ad_arg = sun_dec_deg.to_radians().tan() * lat_deg.to_radians().tan();
    let ascensional_diff = if ad_arg >= 1.0 {
        90.0 - SUNSHINE_NEAR_ZERO_DEG
    } else if ad_arg <= -1.0 {
        -90.0 + SUNSHINE_NEAR_ZERO_DEG
    } else {
        ad_arg.asin().to_degrees()
    };
    let nocturnal_semi_arc = 90.0 - ascensional_diff;
    let diurnal_semi_arc = 90.0 + ascensional_diff;

    let lat_rad = lat_deg.to_radians();
    let (sin_lat, cos_lat) = sin_cos(lat_rad);
    let sun_dec_rad = sun_dec_deg.to_radians();
    let cos_sun_dec = sun_dec_rad.cos();
    let tan_sun_dec = sun_dec_rad.tan();
    let eps_rad = eps_deg.to_radians();
    let (sin_obliquity, cos_obliquity) = sin_cos(eps_rad);

    let nsa_third = nocturnal_semi_arc / 3.0;
    const NSA_OFFSETS: [(usize, f64); 4] = [(2, -2.0), (3, -1.0), (5, 1.0), (6, 2.0)];
    for (house_index, multiple) in NSA_OFFSETS {
        cusps[house_index] = sunshine_arc_to_ecliptic(
            multiple * nsa_third,
            false,
            armc_deg,
            lat_deg,
            sun_dec_deg,
            sin_lat,
            cos_lat,
            cos_sun_dec,
            tan_sun_dec,
            sin_obliquity,
            cos_obliquity,
        );
    }
    let dsa_third = diurnal_semi_arc / 3.0;
    const DSA_OFFSETS: [(usize, f64); 4] = [(8, -2.0), (9, -1.0), (11, 1.0), (12, 2.0)];
    for (house_index, multiple) in DSA_OFFSETS {
        cusps[house_index] = sunshine_arc_to_ecliptic(
            multiple * dsa_third,
            true,
            armc_deg,
            lat_deg,
            sun_dec_deg,
            sin_lat,
            cos_lat,
            cos_sun_dec,
            tan_sun_dec,
            sin_obliquity,
            cos_obliquity,
        );
    }

    if mc_under_horizon {
        for i in [2_usize, 3, 5, 6, 8, 9, 11, 12] {
            cusps[i] = wrap360(cusps[i] + 180.0);
        }
    }
    strip_index0(cusps)
}

/// Makransky ecliptic cut: the
/// compact Cartesian/algebraic form of the published *Solar House Cusp
/// Algorithm* steps 6–9 (Makransky 1990, pp. 147–148),
/// `tan λ = sin W / (cos ε·cos W − sin ε·tan p)` with `atan2` quadrant
/// resolution. `W` is the meridian-relative right-ascension blend in
/// degrees, `p` the pole height in degrees, `eps_rad` the obliquity in
/// radians. The sign of `p` selects the cusp-number family
/// (`eastern_family`: cusps 11/12/2/3 positive, 5/6/8/9 negative); the
/// magnitude (not the southern-hemisphere algebraic sign) is
/// intentional. Returns degrees in `[0, 360)`; scalar on both sides.
pub(super) fn makransky_ecliptic_cut(
    w_deg: f64,
    pole_deg: f64,
    eps_rad: f64,
    eastern_family: bool,
) -> f64 {
    let effective_pole = if eastern_family {
        pole_deg.abs()
    } else {
        -pole_deg.abs()
    };
    let w_rad = w_deg.to_radians();
    let pole_rad = effective_pole.to_radians();
    let (sin_w, cos_w) = sin_cos(w_rad);
    let (sin_eps, cos_eps) = sin_cos(eps_rad);
    let longitude = sin_w.atan2(cos_eps * cos_w - sin_eps * pole_rad.tan());
    wrap360(longitude.to_degrees())
}

/// One non-angular Sunshine cusp from Makransky's construction
/// RA offset `X` from the
/// relevant meridian (night-arc offsets from RAIC, day-arc from RAMC),
/// meridian-distance continuation across the east/west hour circle
/// (`180 − |X|`, upper/lower swap, `|A + F − 90°|`, `W = RA ∓ Q` sign
/// reversal), the `MD = 90°` limiting formula, and the `|X| > 90°`
/// zenith-distance reduction — the standard spherical-coordinate
/// continuation of the printed acute-triangle algorithm (Makransky
/// 1988 pp. 36–39, 1990 pp. 145–148, fig. IV-4).
///
/// `index` is one of the eight intermediate cusps 2, 3, 5, 6, 8, 9, 11,
/// 12 (nocturnal 2/3/5/6 take night-arc offsets from RAIC, day arcs
/// from RAMC); `offset_deg`, `armc_deg`, `lat_deg`, `sun_dec_deg`
/// are degrees, `eps_rad` radians. `armc_deg` is the working
/// upper-meridian right ascension (RAIC already substituted for
/// southern latitudes by the caller). Returns
/// degrees in `[0, 360)`; scalar on both sides.
pub(super) fn makransky_intermediate_cusp(
    index: usize,
    offset_deg: f64,
    armc_deg: f64,
    lat_deg: f64,
    sun_dec_deg: f64,
    eps_rad: f64,
) -> f64 {
    let nocturnal = matches!(index, 2 | 3 | 5 | 6);
    let eastern_family = matches!(index, 2 | 3 | 11 | 12);

    let ra_deg = wrap360(armc_deg + if nocturnal { 180.0 } else { 0.0 } + offset_deg);
    let raw_md = offset_deg.abs();
    let beyond_east_west = raw_md > 90.0;
    let md_deg = if beyond_east_west {
        180.0 - raw_md
    } else {
        raw_md
    };

    let is_upper_distance = (!nocturnal) != beyond_east_west;
    let abs_lat_deg = lat_deg.abs();
    let lat_rad = lat_deg.to_radians();
    let abs_lat_rad = abs_lat_deg.to_radians();
    let dec_rad = sun_dec_deg.to_radians();
    let md_rad = md_deg.to_radians();
    let (sin_lat, cos_lat) = sin_cos(lat_rad);
    let sin_abs_lat = abs_lat_rad.sin();
    let (sin_md, cos_md) = sin_cos(md_rad);

    let zenith_distance = if (md_deg - 90.0).abs() < 1e-12 {
        90.0 - (sin_abs_lat * dec_rad.tan()).atan().to_degrees()
    } else {
        let angle_a = (cos_lat * md_rad.tan()).atan().to_degrees();
        let angle_b = (abs_lat_rad.tan() * cos_md).atan().to_degrees();
        let same_declination_hemisphere = (sun_dec_deg >= 0.0) == (lat_deg >= 0.0);
        let c_sign = if same_declination_hemisphere == is_upper_distance {
            -1.0
        } else {
            1.0
        };
        let triangle_c = angle_b + c_sign * sun_dec_deg.abs();
        let angle_f = (sin_abs_lat * sin_md * triangle_c.to_radians().tan())
            .atan()
            .to_degrees();
        let mut zd = angle_a + angle_f;
        if beyond_east_west {
            zd = (zd - 90.0).abs();
        }
        zd
    };

    let pole_argument = sin_lat * zenith_distance.to_radians().sin();
    let pole_deg = clamp_unit(pole_argument).asin().to_degrees();
    let q_argument = dec_rad.tan() * pole_deg.to_radians().tan();
    let q_deg = clamp_unit(q_argument).asin().to_degrees();

    let mut q_sign = if eastern_family { -1.0 } else { 1.0 };
    if beyond_east_west {
        q_sign = -q_sign;
    }
    let w_deg = wrap360(ra_deg + q_sign * q_deg);
    makransky_ecliptic_cut(w_deg, pole_deg, eps_rad, eastern_family)
}

/// Sunshine/Makransky houses: the
/// Sun's diurnal arc in six equal parts and its nocturnal arc in six,
/// the division points on the Sun's parallel carried to the ecliptic
/// along house circles through the horizon north/south points
/// (Makransky, *Solar House System*, pp. 147–148; diurnal/nocturnal
/// angles 1/4/7/10). The four angles are the horizon/meridian circles;
/// opposite intermediate cusps are generally not opposite in longitude.
///
/// Arguments are degrees (active ARMC, latitude, true obliquity,
/// precomputed Ascendant/MC, `sun_dec_deg` explicit — the release
/// `ascmc9` slot; entry owns its acquisition). Takes the passed angles
/// as-is. Returns `[f64; 12]` houses 1–12.
/// Southern latitudes use RAIC in place of RAMC with a 180° output
/// rotation (Makransky p. 148; exact coordinate transformation, the
/// supplied angular cusps already hemisphere-correct). Exact zero
/// `sun_dec`/`lat` evaluate through the same equations (no sign-bit
/// branches), so the equinox carries no 180° jump.
///
/// Returns [`Error::PolarCircle`] with `house_system == 'i'` when the
/// Sun is circumpolar (`|tan δ☉·tan φ| >= 1`, gate) — the
/// release `PolarCircleError`. The threshold is the exact circumpolar
/// latitude bound `90 − |sun_dec|`
/// (`|φ| >= 90 − |δ| ⟺ |tan δ·tan φ| >= 1` for `|φ|, |δ| < 90°`);
/// `sun_dec` is the release `ascmc9` slot (entry acquisition stays with
/// the ARMC-entry layer). Its exact consumed word also belongs to the owned
/// constructor presentation; the typed primary retains the SDK threshold.
pub fn houses_sunshine_makransky(
    armc_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
    asc_deg: f64,
    mc_deg: f64,
    sun_dec_deg: f64,
) -> Result<[f64; 12]> {
    let tan_product = sun_dec_deg.to_radians().tan() * lat_deg.to_radians().tan();
    if tan_product.abs() >= 1.0 {
        return Err(Error::PolarCircle {
            latitude: lat_deg,
            threshold: 90.0 - sun_dec_deg.abs(),
            obliquity: eps_deg,
            house_system: 'i',
        });
    }
    let ascensional_difference = tan_product.asin().to_degrees();
    let nocturnal_semiarc = 90.0 - ascensional_difference;
    let diurnal_semiarc = 90.0 + ascensional_difference;
    let night_third = nocturnal_semiarc / 3.0;
    let day_third = diurnal_semiarc / 3.0;

    let southern = lat_deg < 0.0;
    let working_armc = if southern {
        wrap360(armc_deg + 180.0)
    } else {
        armc_deg
    };
    let output_rotation = if southern { 180.0 } else { 0.0 };
    let eps_rad = eps_deg.to_radians();

    let mut cusps = [0.0_f64; 13];
    cusps[1] = asc_deg;
    cusps[7] = wrap360(asc_deg + 180.0);
    cusps[10] = mc_deg;
    cusps[4] = wrap360(mc_deg + 180.0);
    const OFFSETS: [(usize, f64, bool); 8] = [
        (2, -2.0, false),
        (3, -1.0, false),
        (5, 1.0, false),
        (6, 2.0, false),
        (8, -2.0, true),
        (9, -1.0, true),
        (11, 1.0, true),
        (12, 2.0, true),
    ];
    for (cusp, multiple, diurnal) in OFFSETS {
        let offset = multiple * if diurnal { day_third } else { night_third };
        let longitude =
            makransky_intermediate_cusp(cusp, offset, working_armc, lat_deg, sun_dec_deg, eps_rad);
        cusps[cusp] = wrap360(longitude + output_rotation);
    }
    Ok(strip_index0(cusps))
}

/// Topocentric/specialized dispatch over the topocentric/specialized-layer
/// systems with the release `houses`/`houses_armc` branch order for these codes.
///
/// `hsys` is folded here ([`fold_hsys_case`]), like the release folds
/// before dispatching. `sun_dec` serves only `I`/`i`; every other kernel ignores it. The R/C/T/J
/// MC-flip prologue selecting `armc_active`/`calc_lat`, the APC/I/i
/// `ascmc[1]` re-sync and the full Ascendant/MC acquisition stay with the
/// ARMC-entry scope — the kernels take the passed angles as-is.
///
/// Every other selector — the simple-ecliptic-layer systems, the
/// quadrant-layer systems (which keep their own [`houses_quadrant`]
/// dispatcher, including the unknown-selector Placidus default and the
/// served Koch `K` arm) and Gauquelin (`G`/`g`, Gauquelin-sector layer) — fails with
/// [`Error::HouseSystemUnserved`],
/// never with a silent substitution. None of the topocentric/specialized-layer
/// systems takes the polar-circle refusal; the Makransky circumpolar-Sun refusal
/// surfaces from [`houses_sunshine_makransky`].
pub fn houses_specialized(
    armc_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
    asc_deg: f64,
    mc_deg: f64,
    sun_dec_deg: f64,
    hsys: char,
) -> Result<[f64; 12]> {
    let code = fold_hsys_case(hsys);
    match code {
        'T' => Ok(houses_polich_page(
            armc_deg, lat_deg, eps_deg, asc_deg, mc_deg,
        )),
        'H' => Ok(houses_horizontal(
            armc_deg, lat_deg, eps_deg, asc_deg, mc_deg,
        )),
        'Y' => Ok(houses_apc(armc_deg, lat_deg, eps_deg, asc_deg, mc_deg)),
        'F' => Ok(houses_carter(armc_deg, lat_deg, eps_deg, asc_deg, mc_deg)),
        'U' => houses_krusinski(armc_deg, lat_deg, eps_deg, asc_deg, mc_deg),
        'N' => Ok(houses_natural_gradient(
            armc_deg, lat_deg, eps_deg, asc_deg, mc_deg,
        )),
        'I' => Ok(houses_sunshine(
            armc_deg,
            lat_deg,
            eps_deg,
            asc_deg,
            mc_deg,
            sun_dec_deg,
        )),
        'i' => houses_sunshine_makransky(armc_deg, lat_deg, eps_deg, asc_deg, mc_deg, sun_dec_deg),
        'J' => Ok(houses_savard_a(armc_deg, lat_deg, eps_deg, asc_deg, mc_deg)),
        other => Err(Error::HouseSystemUnserved {
            house_system: other,
        }),
    }
}

/// One intermediate Gauquelin sector cusp via iterative pole-height
/// refinement (Gauquelin 1955,
/// Meeus 2nd ed., Ch. 13): horizon-fraction `sector_offset/9` of the
/// ascensional-difference arc projected onto the ecliptic, with the
/// spherical-trigonometry core shared by both quadrant loops.
///
/// `sector_offset` is the integer horizon-fraction index (`1–8`) of the
/// quadrant arc; `base_ramc_deg` is the quadrant base in degrees
/// (`armc` for the fourth/second-quadrant loop, `armc + 180` for the
/// first/third) and `ramc_sign` its direction (`+1.0`/`−1.0`);
/// `ascensional_diff_deg`, `eps_deg`, `lat_deg` are degrees and the
/// obliquity/latitude trig legs pre-computed so the sixteen sector
/// calls share them, exactly as the release builds one `_cusp_kwargs`
/// dict. Returns one ecliptic longitude in degrees `[0, 360)`; scalar
/// on both sides (no sector-layout line applies here).
#[allow(clippy::too_many_arguments)]
pub(super) fn gauquelin_cusp_for_sector(
    sector_offset: u32,
    base_ramc_deg: f64,
    ramc_sign: f64,
    ascensional_diff_deg: f64,
    tan_obliquity: f64,
    sin_obliquity: f64,
    tan_lat: f64,
    eps_deg: f64,
    lat_deg: f64,
) -> f64 {
    // Intermediate RAMC for this sector.
    let intermediate_ramc =
        wrap360(base_ramc_deg + ramc_sign * (90.0 / 9.0) * f64::from(sector_offset));

    // Degenerate obliquity (eps -> 0): the ecliptic coincides with the
    // equator, so every sector cusp lies on the equator at the
    // intermediate RAMC — the same value the equatorial-declination
    // branch below returns. Guards the division by tan(eps).
    if tan_obliquity.abs() < GAUQUELIN_NEAR_ZERO {
        return intermediate_ramc;
    }

    // Initial pole height: latitude of the great-circle pole for this
    // sector fraction.
    let arc_fraction_deg = ascensional_diff_deg * f64::from(sector_offset) / 9.0;
    let initial_pole_height = (arc_fraction_deg.to_radians().sin() / tan_obliquity)
        .atan()
        .to_degrees();

    // Seed estimate: ecliptic point rising at the intermediate RAMC with
    // the initial pole.
    let mut cusp = calc_ascendant(intermediate_ramc, eps_deg, lat_deg, initial_pole_height);

    // Declination of the seed cusp point.
    let cusp_dec = (sin_obliquity * cusp.to_radians().sin()).asin();
    let tan_cusp_dec = cusp_dec.tan();

    // Degenerate case: cusp on the equator, no pole-height correction
    // possible.
    if tan_cusp_dec.abs() < GAUQUELIN_NEAR_ZERO {
        return intermediate_ramc;
    }

    // First pole-height refinement from the seed (
    // `degrees(atan(sin(radians(...)) / tan_cusp_dec))` association; a
    // plain quotient-atan, not atan2, so negative denominators keep the
    // release branch).
    let pole_arg = clamp_unit(tan_lat * tan_cusp_dec);
    let pole_height = ((pole_arg.asin().to_degrees() * f64::from(sector_offset) / 9.0)
        .to_radians()
        .sin()
        / tan_cusp_dec)
        .atan()
        .to_degrees();
    cusp = calc_ascendant(intermediate_ramc, eps_deg, lat_deg, pole_height);

    // Iterative refinement until convergence.
    let mut prev_cusp = 0.0;
    for iteration in 1..=GAUQUELIN_MAX_ITERATIONS {
        let updated_dec = (sin_obliquity * cusp.to_radians().sin()).asin();
        let tan_cusp_dec_updated = updated_dec.tan();

        if tan_cusp_dec_updated.abs() < GAUQUELIN_NEAR_ZERO {
            cusp = intermediate_ramc;
            break;
        }

        let pole_arg = clamp_unit(tan_lat * tan_cusp_dec_updated);
        let pole_height = ((pole_arg.asin().to_degrees() * f64::from(sector_offset) / 9.0)
            .to_radians()
            .sin()
            / tan_cusp_dec_updated)
            .atan()
            .to_degrees();
        cusp = calc_ascendant(intermediate_ramc, eps_deg, lat_deg, pole_height);

        if iteration > 1 {
            let mut angular_diff = (cusp - prev_cusp).abs();
            if angular_diff > 180.0 {
                angular_diff = 360.0 - angular_diff;
            }
            if angular_diff < GAUQUELIN_CONVERGENCE_DEG {
                break;
            }
        }
        prev_cusp = cusp;
    }

    cusp
}

/// Gauquelin 36-sector houses:
/// 36 sectors of 10° — the Placidus-inspired diurnal/nocturnal
/// semi-arc ninths projected onto the ecliptic via pole-height
/// refinement (Gauquelin 1955; Meeus 2nd ed., Ch. 13).
///
/// With `a = arcsin(tan φ·tan ε)`: sectors 2–9 take offsets
/// `+80°…+10°` (`10 − ih`) from the active ARMC, sectors 29–36 take
/// offsets `−10°…−80°` (`ih − 28`) from `ARMC + 180°`, each refined
/// through `gauquelin_cusp_for_sector`; opposites at `+180°`;
/// cardinals pinned (sector 1 = supplied Ascendant, 10 = supplied MC,
/// 19 = Desc, 28 = IC). Arguments are degrees (active ARMC — the entry
/// owns the polar-circle refusal and Ascendant/MC acquisition; this
/// kernel takes them as given — latitude, true obliquity, supplied
/// Ascendant/MC cardinals).
///
/// Inside the polar-circle fallback band (`|lat| >= 90 − eps`) the
/// sectors are equal 10° divisions counting down from the supplied
/// Ascendant. Returns `[f64; 36]` sectors 1–36 in sector order (the
/// release 37-list with index 0 unused, rendered in Rust order); the
/// `g` 12-sector prefix is [`houses_gauquelin_12`].
pub fn houses_gauquelin(
    armc_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
    asc_deg: f64,
    mc_deg: f64,
) -> [f64; 36] {
    // Index 0 unused, sectors 1–36.
    let mut cusps = [0.0_f64; 37];

    // Polar-circle fallback: equal 10° divisions from the Ascendant.
    if lat_deg.abs() >= 90.0 - eps_deg {
        for (i, cusp) in cusps.iter_mut().enumerate().skip(1) {
            *cusp = wrap360(asc_deg - (i as f64 - 1.0) * 10.0);
        }
        return strip_gauquelin_index0(cusps);
    }

    // Precomputed obliquity and latitude trig quantities.
    let eps_rad = eps_deg.to_radians();
    let (sin_obliquity, _) = sin_cos(eps_rad);
    let tan_obliquity = eps_rad.tan();
    let tan_lat = lat_deg.to_radians().tan();

    // Ascensional difference: a = arcsin(tan(φ) · tan(ε)).
    let asc_diff_arg = clamp_unit(tan_lat * tan_obliquity);
    let ascensional_diff = asc_diff_arg.asin().to_degrees();

    // Fourth/second quadrant: sectors 2–9 (RAMC offsets +10° … +80°).
    for ih in 2..10_u32 {
        let sector_offset = 10 - ih;
        cusps[ih as usize] = gauquelin_cusp_for_sector(
            sector_offset,
            armc_deg,
            1.0,
            ascensional_diff,
            tan_obliquity,
            sin_obliquity,
            tan_lat,
            eps_deg,
            lat_deg,
        );
        cusps[(ih + 18) as usize] = wrap360(cusps[ih as usize] + 180.0);
    }

    // First/third quadrant: sectors 29–36 (RAMC offsets −10° … −80°).
    for ih in 29..37_u32 {
        let sector_offset = ih - 28;
        cusps[ih as usize] = gauquelin_cusp_for_sector(
            sector_offset,
            wrap360(armc_deg + 180.0),
            -1.0,
            ascensional_diff,
            tan_obliquity,
            sin_obliquity,
            tan_lat,
            eps_deg,
            lat_deg,
        );
        cusps[(ih - 18) as usize] = wrap360(cusps[ih as usize] + 180.0);
    }

    // Cardinal sectors.
    cusps[1] = asc_deg;
    cusps[10] = mc_deg;
    cusps[19] = wrap360(asc_deg + 180.0);
    cusps[28] = wrap360(mc_deg + 180.0);

    strip_gauquelin_index0(cusps)
}

/// Strip the unused index 0 of the release 37-slot Gauquelin layout to
/// the public 36-sector order (sectors 1–36).
pub(super) fn strip_gauquelin_index0(cusps: [f64; 37]) -> [f64; 36] {
    let mut out = [0.0_f64; 36];
    out.copy_from_slice(&cusps[1..37]);
    out
}

/// Gauquelin `g` return shape (`houses`/`houses_armc` rule):
/// the same 36-sector computation, returning the first 12 sectors in
/// the ordinary 12-cusp tuple.
pub fn houses_gauquelin_12(
    armc_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
    asc_deg: f64,
    mc_deg: f64,
) -> [f64; 12] {
    let full = houses_gauquelin(armc_deg, lat_deg, eps_deg, asc_deg, mc_deg);
    let mut out = [0.0_f64; 12];
    out.copy_from_slice(&full[..12]);
    out
}

/// Compute twelve Koch house cusps from the ascensional-difference construction.
///
/// `armc_deg` is right ascension of the meridian, `lat_deg` is geographic
/// latitude, `eps_deg` is obliquity, and `asc_deg`/`mc_deg` are the previously
/// computed ascendant and midheaven. All inputs and returned ecliptic longitudes
/// are in degrees; array index zero denotes house 1. The construction clamps
/// inverse-trigonometric arguments and sets opposite cusps 180 degrees apart.
///
/// For `abs(lat_deg) >= 90 - eps_deg`, returns Porphyry cusps. This low-level
/// kernel does not validate inputs or return a polar warning; callers must
/// supply finite, mutually consistent angles and a physical latitude/obliquity.
pub fn houses_koch(
    armc_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
    asc_deg: f64,
    mc_deg: f64,
) -> [f64; 12] {
    // Kernel-level polar fallback ( `>=` seam): Porphyry cusps.
    if lat_deg.abs() >= 90.0 - eps_deg {
        let porphyry = houses_porphyry(asc_deg, mc_deg);
        let mut full = [0.0_f64; 13];
        full[1..13].copy_from_slice(&porphyry);
        return strip_index0(full);
    }

    // Precomputed obliquity and latitude trig quantities .
    let sin_obliquity = eps_deg.to_radians().sin();
    let tan_lat = lat_deg.to_radians().tan();
    let cos_lat = lat_deg.to_radians().cos();

    // Koch auxiliary: sin(dec_MC) = sin(MC)·sin(ε)/cos(φ), clamped.
    // The latitude cosine belongs to the Koch/Knappich construction
    // (see Holden 1977), not to Meeus.
    let mc_dec_sin = clamp_unit(mc_deg.to_radians().sin() * sin_obliquity / cos_lat);
    // Always >= 0.
    let mc_dec_cos = (1.0 - mc_dec_sin * mc_dec_sin).sqrt();

    // Auxiliary angle c: elevation of the MC above the prime vertical
    // (`atan2` handles the degenerate `mc_dec_cos == 0` natively).
    let zenith_angle = tan_lat.atan2(mc_dec_cos).to_degrees();

    // Ascensional difference divided by 3 (offsets each intermediate cusp).
    let asc_diff_arg = clamp_unit(zenith_angle.to_radians().sin() * mc_dec_sin);
    let asc_diff_third = asc_diff_arg.asin().to_degrees() / 3.0;

    // Intermediate cusps are ecliptic rising points at ARMC offsets of
    // 30, 60, 120, 150 degrees, adjusted by multiples of the Koch arc.
    let mut cusps = init_cardinal_cusps_13(asc_deg, mc_deg);
    cusps[11] = calc_ascendant(
        armc_deg + 30.0 - 2.0 * asc_diff_third,
        eps_deg,
        lat_deg,
        lat_deg,
    );
    cusps[12] = calc_ascendant(armc_deg + 60.0 - asc_diff_third, eps_deg, lat_deg, lat_deg);
    cusps[2] = calc_ascendant(armc_deg + 120.0 + asc_diff_third, eps_deg, lat_deg, lat_deg);
    cusps[3] = calc_ascendant(
        armc_deg + 150.0 + 2.0 * asc_diff_third,
        eps_deg,
        lat_deg,
        lat_deg,
    );

    // Opposite houses.
    set_opposite_cusps_13(&mut cusps);
    strip_index0(cusps)
}

// ASCMC/auxiliary angles and the `houses_armc` entry.

/// Vertex: the auxiliary western
/// prime-vertical/ecliptic intersection (fateful-encounter use), from
/// `num = −cos(ARMC)`, `den = sin(ARMC)·cos ε − sin ε/tan φ` via
/// `atan2`.
///
/// `armc_deg` is the true (unflipped) ARMC, `eps_deg` the true
/// obliquity, `lat_deg` the geographic latitude — all degrees in,
/// longitude out in `[0, 360)`. At exactly `lat == 0.0` the
/// `1/tan(lat)` singularity takes the deterministic positive-latitude
/// one-sided limit (`1e-12`); a joint `1e-12` numerator/denominator
/// degeneracy (latitude on the solstitial colure crossing, `|lat| ==
/// eps` with ARMC 90/270) pins the western cardinal `armc + 270`. The
/// `atan2` candidate is disambiguated to the western antipode in the
/// equatorial frame (candidate-RA hour angle in `(0, 180)`, else
/// `+180°`). Precision is ~0.001° away from the equator. Scalar on both
/// sides.
pub fn calc_vertex(armc_deg: f64, eps_deg: f64, lat_deg: f64) -> f64 {
    let eps_rad = eps_deg.to_radians();
    let (sin_eps, cos_eps) = sin_cos(eps_rad);
    let lat = if lat_deg == 0.0 { 1e-12 } else { lat_deg };
    let armc_rad = armc_deg.to_radians();
    let lat_rad = lat.to_radians();
    let (sin_armc, cos_armc) = sin_cos(armc_rad);
    let num = -cos_armc;
    let den = sin_armc * cos_eps - sin_eps / lat_rad.tan();
    if num.abs() < 1e-12 && den.abs() < 1e-12 {
        return wrap360(armc_deg + 270.0);
    }
    let mut vtx = wrap360(num.atan2(den).to_degrees());
    let vtx_rad = vtx.to_radians();
    let (sin_vtx, cos_vtx) = sin_cos(vtx_rad);
    let ra_vtx = wrap360((cos_eps * sin_vtx).atan2(cos_vtx).to_degrees());
    let hour_angle = wrap360(armc_deg - ra_vtx);
    if !(0.0 < hour_angle && hour_angle < 180.0) {
        vtx = wrap360(vtx + 180.0);
    }
    vtx
}

/// Medium Coeli from the (possibly flipped) ARMC: the meridian/ecliptic culmination from the sidereal
/// angle, `MC = atan2(tan(ARMC), cos ε)`, with the quadrant correction
/// keeping MC in the same ecliptic half as the ARMC (the `(90, 270]`
/// vs rest rule with the `<=` boundary matching the canonical entry),
/// then `% 360`.
///
/// Arguments are degrees in, longitude out in `[0, 360)`; scalar on
/// both sides.
pub fn mc_from_armc(armc_active_deg: f64, eps_deg: f64) -> f64 {
    let eps_rad = eps_deg.to_radians();
    let mut mc = armc_active_deg
        .to_radians()
        .tan()
        .atan2(eps_rad.cos())
        .to_degrees();
    if mc < 0.0 {
        mc += 360.0;
    }
    if 90.0 < armc_active_deg && armc_active_deg <= 270.0 {
        if mc <= 90.0 || mc > 270.0 {
            mc += 180.0;
        }
    } else if armc_active_deg > 270.0 {
        if mc <= 270.0 {
            mc += 180.0;
        }
    // Collapsed from the release's nested `elif armc <= 90: if mc > 90`
    // (identical condition order and NaN behavior; clippy collapsible-if).
    } else if armc_active_deg <= 90.0 && mc > 90.0 {
        mc += 180.0;
    }
    wrap360(mc)
}

/// Ascendant from the true ARMC (Smart/Meeus rising-point construction): the horizon
/// rising point from `num = cos(ARMC)`,
/// `den = −(sin(ARMC)·cos ε + tan φ·sin ε)` via `atan2`, with
/// eastern-horizon disambiguation via RA/Dec/hour-angle/azimuth against
/// the true ARMC (`+180°` when `0 < H < 180`); the `1e-12` joint
/// degeneracy pins the east point (`armc + 90`); the `1e-10` near-360
/// snap reports `0.0`.
///
/// Arguments are degrees in, longitude out in `[0, 360)`; scalar on
/// both sides.
pub fn asc_from_armc(armc_deg: f64, eps_deg: f64, lat_deg: f64) -> f64 {
    let eps_rad = eps_deg.to_radians();
    let (sin_eps, cos_eps) = sin_cos(eps_rad);
    let armc_rad = armc_deg.to_radians();
    let (sin_armc, cos_armc) = sin_cos(armc_rad);
    let lat_rad = lat_deg.to_radians();
    let (sin_lat, cos_lat) = sin_cos(lat_rad);
    let num = cos_armc;
    let den = -(sin_armc * cos_eps + lat_rad.tan() * sin_eps);
    let mut asc = if num.abs() < 1e-12 && den.abs() < 1e-12 {
        wrap360(armc_deg + 90.0)
    } else {
        let mut asc = wrap360(num.atan2(den).to_degrees());
        let asc_rad = asc.to_radians();
        let (sin_asc, cos_asc) = sin_cos(asc_rad);
        let ra = wrap360((cos_eps * sin_asc).atan2(cos_asc).to_degrees());
        let dec_rad = (sin_eps * sin_asc).asin();
        let h_deg = wrap360(armc_deg - ra + 360.0);
        let h_rad = h_deg.to_radians();
        let (sin_h, cos_h) = sin_cos(h_rad);
        let num_az = sin_h;
        let den_az = sin_lat * cos_h - cos_lat * dec_rad.tan();
        let _az = wrap360(num_az.atan2(den_az).to_degrees() + 180.0);
        if 0.0 < h_deg && h_deg < 180.0 {
            asc = wrap360(asc + 180.0);
        }
        asc
    };
    if asc >= 360.0 || (360.0 - asc) < 1e-10 {
        asc = 0.0;
    }
    asc
}

/// Equatorial Ascendant / East Point: the ecliptic longitude whose right ascension is
/// `ARMC + 90°`, `λ = atan2(sin RA, cos RA·cos ε)` with
/// `RA = (ARMC + 90) mod 360`. Arguments are degrees in, longitude out
/// in `[0, 360)`; scalar on both sides.
pub fn equatorial_ascendant_deg(armc_deg: f64, eps_deg: f64) -> f64 {
    let equ_asc_ra = wrap360(armc_deg + 90.0);
    let equ_asc_ra_rad = equ_asc_ra.to_radians();
    let eps_rad = eps_deg.to_radians();
    let (sin_ra, cos_ra) = sin_cos(equ_asc_ra_rad);
    wrap360((sin_ra).atan2(cos_ra * eps_rad.cos()).to_degrees())
}

/// Co-Ascendant W. Koch: the
/// antipode of the `ARMC − 90°` rising point,
/// `coasc1 = Asc(ARMC − 90°, φ) + 180°`. Arguments are degrees in,
/// longitude out in `[0, 360)`; scalar on both sides.
pub fn coasc_koch_deg(armc_deg: f64, eps_deg: f64, lat_deg: f64) -> f64 {
    let coasc_armc = wrap360(armc_deg - 90.0);
    wrap360(calc_ascendant(coasc_armc, eps_deg, lat_deg, lat_deg) + 180.0)
}

/// Co-Ascendant M. Munkasey: the
/// rising point at `ARMC + 90°` with the complementary pole height
/// (`90 − φ` for `lat > 0`, `−90 − φ` otherwise). At exactly
/// `lat == 0.0` the horizon system `H` takes the southern branch and
/// every other system the northern branch (per-system exact-equator
/// rule, not the short `lat >= 0` comment form). `hsys` is a branch
/// selector only. Arguments are degrees in, longitude out in `[0, 360)`;
/// scalar on both sides.
pub fn coasc_munkasey_deg(armc_deg: f64, eps_deg: f64, lat_deg: f64, hsys: char) -> f64 {
    let coasc2_armc = wrap360(armc_deg + 90.0);
    let coasc2_lat = if lat_deg > 0.0 || (lat_deg == 0.0 && hsys != 'H') {
        90.0 - lat_deg
    } else {
        -90.0 - lat_deg
    };
    calc_ascendant(coasc2_armc, eps_deg, coasc2_lat, coasc2_lat)
}

/// Polar Ascendant M. Munkasey:
/// `polasc = Asc(ARMC − 90°, φ)` — the same rising point as coasc1
/// without the Koch `+180°`. Arguments are degrees in, longitude out in
/// `[0, 360)`; scalar on both sides.
pub fn polar_ascendant_deg(armc_deg: f64, eps_deg: f64, lat_deg: f64) -> f64 {
    let coasc_armc = wrap360(armc_deg - 90.0);
    calc_ascendant(coasc_armc, eps_deg, lat_deg, lat_deg)
}

/// Full 8-slot ASCMC assembly:
/// the reference-API 8-angle vector — `[Asc, MC, ARMC, Vertex, EquAsc,
/// CoAscKoch, CoAscMunkasey, PolarAsc]`: the rising point, the meridian
/// culmination, the sidereal-time angle, the western prime-vertical
/// intersection, the east point, the Koch +180° antipode, the Munkasey
/// complementary-height point, the unflipped Koch-base point — with the
/// MC computed from the flipped `armc_active` and every other angle
/// from the true `armc_deg`.
///
/// All angles are degrees in; `hsys` only steers the coasc2 equatorial
/// branch. Returns `[f64; 8]` scalar degrees (the assembly, not cusps).
pub fn compute_ascmc(
    armc_deg: f64,
    armc_active_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
    hsys: char,
) -> [f64; 8] {
    let mut ascmc = [0.0_f64; N_ASCMC];
    ascmc[ASCMC_ASC] = asc_from_armc(armc_deg, eps_deg, lat_deg);
    ascmc[ASCMC_MC] = mc_from_armc(armc_active_deg, eps_deg);
    ascmc[ASCMC_ARMC] = armc_deg;
    ascmc[ASCMC_VERTEX] = calc_vertex(armc_deg, eps_deg, lat_deg);
    ascmc[ASCMC_EQUASC] = equatorial_ascendant_deg(armc_deg, eps_deg);
    ascmc[ASCMC_COASC_KOCH] = coasc_koch_deg(armc_deg, eps_deg, lat_deg);
    ascmc[ASCMC_COASC_MUNKASEY] = coasc_munkasey_deg(armc_deg, eps_deg, lat_deg, hsys);
    ascmc[ASCMC_POLAR_ASC] = polar_ascendant_deg(armc_deg, eps_deg, lat_deg);
    ascmc
}

/// Whether the MC sits below the horizon for the R/C/T/J flip and the
/// APC re-sync (`houses`/`houses_armc` expression):
/// `sin_alt = sin φ·sin δ_MC + cos φ·cos δ_MC < 0` with
/// `δ_MC = atan(sin(ARMC)·tan ε)` evaluated at the true ARMC.
pub(super) fn mc_below_horizon(armc_deg: f64, eps_deg: f64, lat_deg: f64) -> bool {
    let armc_rad = armc_deg.to_radians();
    let eps_rad = eps_deg.to_radians();
    let mc_dec_rad = (armc_rad.sin() * eps_rad.tan()).atan();
    let lat_rad = lat_deg.to_radians();
    let (sin_lat, cos_lat) = sin_cos(lat_rad);
    let (sin_mcdec, cos_mcdec) = sin_cos(mc_dec_rad);
    sin_lat * sin_mcdec + cos_lat * cos_mcdec < 0.0
}

/// Flipped ARMC for the R/C/T/J MC-below-horizon rule: keep the MC above the horizon for the
/// R/C/T/J polar-circle constructions — `δ_MC = atan(sin(ARMC)·tan ε)`,
/// `sin_alt = sin φ·sin δ_MC + cos φ·cos δ_MC`, and `sin_alt < 0`
/// (strict: no equality flip) reports `armc + 180`, otherwise the true
/// ARMC. Arguments are degrees in (`hsys` gates the flip to R/C/T/J),
/// flipped-or-true ARMC out; scalar on both sides.
pub fn resolve_armc_active(armc_deg: f64, eps_deg: f64, lat_deg: f64, hsys: char) -> f64 {
    if matches!(hsys, 'R' | 'C' | 'T' | 'J') && mc_below_horizon(armc_deg, eps_deg, lat_deg) {
        wrap360(armc_deg + 180.0)
    } else {
        armc_deg
    }
}

/// Latitude for the intermediate-cusp computation (same order as
/// `houses_armc`): negated when the MC flip moved `armc_active` off the
/// true ARMC (exact float inequality, like the release `!=`).
pub fn calc_lat_for_entry(armc_active_deg: f64, armc_deg: f64, lat_deg: f64) -> f64 {
    if armc_active_deg != armc_deg {
        -lat_deg
    } else {
        lat_deg
    }
}

/// Low-precision apparent solar declination (Meeus 1998, ch. 25): the Sunshine `I`/`i`
/// fallback when the ephemeris does not cover the date — avoiding the
/// severe high-latitude distortion of a `0.0` substitute — via TT
/// centuries `t = (tjdut + ΔT − 2451545)/36525`, mean longitude `l0`,
/// anomaly `m`, equation-of-center `c`, node `omega`, apparent longitude
/// `lon_app`, nutation-corrected `eps`, and
/// `δ = asin(sin eps·sin lon_app)`.
///
/// `tjdut` is a Julian-day UT epoch and `deltat_days` the Delta-T value
/// in days; returns declination degrees.
/// Accuracy is ~0.01°, ample for cusp geometry. Scalar on both sides.
pub fn sun_declination_analytic(tjdut: f64, deltat_days: f64) -> f64 {
    let t = (tjdut + deltat_days - 2451545.0) / 36525.0;
    let l0 = wrap360(280.46646 + 36000.76983 * t + 0.0003032 * t * t);
    let m = wrap360(357.52911 + 35999.05029 * t - 0.0001537 * t * t).to_radians();
    let c = (1.914602 - 0.004817 * t - 0.000014 * t * t) * m.sin()
        + (0.019993 - 0.000101 * t) * (2.0 * m).sin()
        + 0.000289 * (3.0 * m).sin();
    let omega = (125.04 - 1934.136 * t).to_radians();
    let (sin_omega, cos_omega) = sin_cos(omega);
    let lon_app = (l0 + c - 0.00569 - 0.00478 * sin_omega).to_radians();
    let eps = (23.4392911111 - 0.0130041667 * t - 1.639e-7 * t * t
        + 5.036e-7 * t * t * t
        + 0.00256 * cos_omega)
        .to_radians();
    (eps.sin() * lon_app.sin()).asin().to_degrees()
}

/// Latitude validation for the ARMC-direct entry: finite
/// latitude in `[-90, 90]` — NaN explicitly refused, since a bare range
/// test lets it through — mapped to [`Error::InvalidObserver`] with
/// `field == "latitude"` (the `CoordinateError` → `InvalidObserver`
/// contract in [`crate::errors`], message prose binding-side). Scalar
/// gate on both sides.
pub fn validate_houses_armc_latitude(lat_deg: f64) -> Result<()> {
    if !lat_deg.is_finite() || lat_deg < -90.0 || lat_deg > 90.0 {
        return Err(Error::InvalidObserver {
            field: "latitude",
            value: lat_deg,
        });
    }
    Ok(())
}

/// ARMC-direct house cusps and ASCMC: cusps plus
/// the ASCMC vector from already-acquired angles, with no JD/ephemeris
/// read.
///
/// Takes the pre-computed ARMC (degrees, normalized `armc % 360`),
/// geographic latitude (validated to `[-90, 90]`), true obliquity
/// (degrees), the folded house-system selector and — for Sunshine `I`/`i`
/// only — the Sun declination in the release `ascmc9` slot (ignored by
/// every other system). The entry prologue owns latitude validation,
/// the R/C/T/J MC-below-horizon flip (`armc_active`/`calc_lat`), the
/// polar-circle refusal for raising systems (raising set per
/// [`polar_raise_applies`], circle per [`is_polar_circle`]), the APC
/// `Y` and `I`/`i` `ascmc[1]` re-syncs (houses 10 and 4 follow the
/// flipped MC), `degnorm` shaping and the Placidus default for
/// truly unknown selectors; the kernels take the active ARMC (plain
/// latitude for B/M/X/H/Y/F/U/N/I/i), the supplied Ascendant/MC
/// cardinals, and for `I`/`i` the explicit `sun_dec` as given. Dispatch
/// serves P/K/R/C/E/A/W/O/B/T/M/X/V/H/F/U/N/S/L/Q/D/I/i/J/Y/G/g in
/// release branch order (beyond the entry docstring's list, which ends
/// at `N`). Returns the cusp longitudes (12 for ordinary systems —
/// the first 12 sectors for lowercase `g` — 36 sectors for uppercase
/// `G`) with the 8 ASCMC angles: Rust `Vec` + `[f64; 8]`, every
/// word `degnorm`-shaped like the release entry.
pub fn houses_armc(
    armc: f64,
    lat_deg: f64,
    eps_deg: f64,
    hsys: char,
    ascmc9_deg: f64,
) -> Result<(Vec<f64>, [f64; N_ASCMC])> {
    houses_armc_buffer(armc, lat_deg, eps_deg, hsys, ascmc9_deg)
        .map(|(cusps, angles)| (cusps.into_vec(), angles))
}

/// Identical ARMC geometry with stack-owned cusp storage for typed SDK calls.
pub(crate) fn houses_armc_buffer(
    armc: f64,
    lat_deg: f64,
    eps_deg: f64,
    hsys: char,
    ascmc9_deg: f64,
) -> Result<(CuspBuffer, [f64; N_ASCMC])> {
    validate_houses_armc_latitude(lat_deg)?;
    let armc_deg = armc.rem_euclid(360.0);
    houses_armc_geometry(armc_deg, lat_deg, eps_deg, hsys, ascmc9_deg)
}

/// Calculate cusps and eight auxiliary angles from an already normalized
/// ARMC in degrees, without applying another remainder operation.
///
/// `armc_deg` is the acquired result of a single remainder by 360 degrees;
/// rounding can leave that result equal to 360. Latitude is checked against
/// the same finite `[-90, 90]` domain as [`houses_armc`]. Obliquity and the
/// Sunshine-only Sun declination are degrees. Output shape, polar refusals
/// and construction errors are identical to [`houses_armc`].
pub fn houses_armc_from_normalized(
    armc_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
    hsys: char,
    ascmc9_deg: f64,
) -> Result<(Vec<f64>, [f64; N_ASCMC])> {
    validate_houses_armc_latitude(lat_deg)?;
    houses_armc_geometry(armc_deg, lat_deg, eps_deg, hsys, ascmc9_deg)
        .map(|(cusps, angles)| (cusps.into_vec(), angles))
}

/// Shared geometry after latitude admission and the caller's one ARMC
/// normalization. No additional remainder is applied at this entry stage.
fn houses_armc_geometry(
    armc_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
    hsys: char,
    ascmc9_deg: f64,
) -> Result<(CuspBuffer, [f64; N_ASCMC])> {
    let key = [
        armc_deg.to_bits(),
        lat_deg.to_bits(),
        eps_deg.to_bits(),
        u64::from(hsys),
        ascmc9_deg.to_bits(),
    ];
    if let Some(geometry) = crate::pure_cache::houses(key) {
        return Ok((CuspBuffer::from(geometry.cusps), geometry.angles));
    }
    let (cusps, angles) = houses_armc_uncached(armc_deg, lat_deg, eps_deg, hsys, ascmc9_deg)?;
    if let Ok(twelve) = cusps.as_slice().try_into() {
        crate::pure_cache::retain_houses(
            key,
            crate::pure_cache::HouseGeometry {
                cusps: twelve,
                angles,
            },
        );
    }
    Ok((cusps, angles))
}

fn houses_armc_uncached(
    armc_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
    hsys: char,
    ascmc9_deg: f64,
) -> Result<(CuspBuffer, [f64; N_ASCMC])> {
    let code = fold_hsys_case(hsys);
    let armc_active = resolve_armc_active(armc_deg, eps_deg, lat_deg, code);
    let mut ascmc = compute_ascmc(armc_deg, armc_active, lat_deg, eps_deg, code);
    let asc = ascmc[ASCMC_ASC];
    let mc = ascmc[ASCMC_MC];
    let calc_lat = calc_lat_for_entry(armc_active, armc_deg, lat_deg);
    if polar_raise_applies(code) && is_polar_circle(lat_deg, eps_deg) {
        return Err(Error::PolarCircle {
            latitude: lat_deg,
            threshold: polar_latitude_threshold_deg(eps_deg),
            obliquity: eps_deg,
            house_system: code,
        });
    }
    // Full dispatch in the release `houses_armc` branch order,
    // including the Placidus default for truly unknown selectors.
    let cusps12: Option<[f64; N_CUSPS]> = match code {
        'P' => Some(houses_placidus(armc_active, lat_deg, eps_deg, asc, mc)),
        'K' => Some(houses_koch(armc_active, lat_deg, eps_deg, asc, mc)),
        'R' => Some(houses_regiomontanus(
            armc_active,
            calc_lat,
            eps_deg,
            asc,
            mc,
        )),
        'C' => Some(houses_campanus(armc_active, calc_lat, eps_deg, asc, mc)),
        'E' | 'A' => Some(houses_equal(asc)),
        'W' => Some(houses_whole_sign(asc)),
        'O' => Some(houses_porphyry(asc, mc)),
        'B' => Some(houses_alcabitius(armc_active, lat_deg, eps_deg, asc, mc)),
        'T' => Some(houses_polich_page(armc_active, calc_lat, eps_deg, asc, mc)),
        'M' => Some(houses_morinus(armc_active, lat_deg, eps_deg, asc, mc)),
        'X' => Some(houses_meridian(armc_active, lat_deg, eps_deg, asc, mc)),
        'V' => Some(houses_vehlow(asc)),
        'H' => Some(houses_horizontal(armc_active, lat_deg, eps_deg, asc, mc)),
        'F' => Some(houses_carter(armc_active, lat_deg, eps_deg, asc, mc)),
        'U' => Some(houses_krusinski(armc_active, lat_deg, eps_deg, asc, mc)?),
        'N' => Some(houses_natural_gradient(
            armc_active,
            lat_deg,
            eps_deg,
            asc,
            mc,
        )),
        'S' => Some(houses_sripati(asc, mc)),
        'L' => Some(houses_pullen_sd(asc, mc)),
        'Q' => Some(houses_pullen_sr(asc, mc)),
        'D' => Some(houses_equal_mc(asc, mc)),
        'I' => {
            let cusps = houses_sunshine(armc_active, lat_deg, eps_deg, asc, mc, ascmc9_deg);
            ascmc[ASCMC_MC] = cusps[9];
            Some(cusps)
        }
        'i' => {
            let cusps =
                houses_sunshine_makransky(armc_active, lat_deg, eps_deg, asc, mc, ascmc9_deg)?;
            ascmc[ASCMC_MC] = cusps[9];
            Some(cusps)
        }
        'J' => Some(houses_savard_a(armc_active, calc_lat, eps_deg, asc, mc)),
        'Y' => {
            let mut cusps = houses_apc(armc_active, lat_deg, eps_deg, asc, mc);
            if mc_below_horizon(armc_deg, eps_deg, lat_deg) {
                ascmc[ASCMC_MC] = wrap360(ascmc[ASCMC_MC] + 180.0);
                cusps[9] = ascmc[ASCMC_MC];
                cusps[3] = wrap360(ascmc[ASCMC_MC] + 180.0);
            }
            Some(cusps)
        }
        'G' | 'g' => None,
        _ => Some(houses_placidus(armc_active, lat_deg, eps_deg, asc, mc)),
    };
    let cusps: CuspBuffer = match cusps12 {
        Some(cusps) => cusps.iter().map(|&c| degnorm_deg(c)).collect(),
        None => {
            let sectors = houses_gauquelin(armc_active, lat_deg, eps_deg, asc, mc);
            if code == 'G' {
                sectors.iter().map(|&c| degnorm_deg(c)).collect()
            } else {
                sectors[..12].iter().map(|&c| degnorm_deg(c)).collect()
            }
        }
    };
    let ascmc = ascmc.map(degnorm_deg);
    Ok((cusps, ascmc))
}

/// Same speed-stencil wrap as `angular_diff_local` / `_rate` /
/// `_rate_richardson` fold): single fold of the raw difference with
/// strict `> 180` / `< −180` guards.
///
/// This is deliberately not [`angular_diff`]: the modulo form maps an
/// exact `−180` difference to `+180`, while the speed fold keeps
/// `−180`. Both agree everywhere else and both feed the same 90°
/// guard, but the fold is kept so the −180° edge carries the
/// release sign.
pub(super) fn wrap_speed_diff_deg(diff: f64) -> f64 {
    if diff > 180.0 {
        diff - 360.0
    } else if diff < -180.0 {
        diff + 360.0
    } else {
        diff
    }
}

/// Centered-difference ARMC rate in deg/day (Meeus ch. 12 sidereal-rate construction):
/// `Δ/(2·d_armc)·_SIDEREAL_RATE` with the divide-then-multiply
/// association, where `d_armc = _SIDEREAL_RATE/86400` is one sidereal
/// second — uniform for every house system — and `_SIDEREAL_RATE =
/// 360.98564736629` ARMC degrees per mean solar day (one extra rotation
/// vs the Sun).
///
/// Wrapped change `> 90°` reports `0.0` (branch discontinuity);
/// sign-pinned `W`/`N` wheels stepping exactly 30° when the Ascendant
/// changes sign report `0.0` on the cusp slots only (`sign_pinned_cusp`
/// selects the cusp rule — the release passes `True` for cusps, default
/// `False` for ASCMC — and the 30° arm additionally requires a
/// sign-pinned system; near-polar genuine ~30° Ascendant sweeps survive
/// on the angle slots, so the rule never touches ASCMC). The selector
/// is folded here (like the release closure folding
/// `chr(_hsys_code(hsys))` on every call), so raw and folded forms
/// behave identically. Arguments are degrees in (`after`/`before`,
/// folded `hsys`, cusp/angle flag), rate in deg/day out (mean solar
/// days via the sidereal rate); scalar on both sides.
pub fn armc_fd_speed_deg_per_day(
    after: f64,
    before: f64,
    sign_pinned_cusp: bool,
    hsys: char,
) -> f64 {
    let diff = wrap_speed_diff_deg(after - before);
    if diff.abs() > SPEED_DISCONTINUITY_DEG {
        return 0.0;
    }
    if sign_pinned_cusp
        && is_sign_pinned(fold_hsys_case(hsys))
        && (diff.abs() - SIGN_PINNED_STEP_DEG).abs() < SIGN_PINNED_STEP_TOL_DEG
    {
        return 0.0;
    }
    diff / (2.0 * ARMC_SPEED_STEP_DEG) * SIDEREAL_RATE_DEG_PER_DAY
}

/// ARMC-direct cusps + ASCMC + their velocities: three full [`houses_armc`] evaluations — center
/// plus ARMC `± d_armc` with the uniform one-sidereal-second stencil
/// for every house system — then the centered-difference rate with
/// `sign_pinned_cusp=True` for cusps, default `False` for ASCMC.
///
/// This is an ARMC-only derivative: obliquity and the supplied Sun declination
/// remain fixed while ARMC advances at [`SIDEREAL_RATE_DEG_PER_DAY`]. It omits
/// their physical time derivatives and any time-dependent ayanamsha; callers
/// with an epoch should prefer [`houses_ex2_via`] for full sampled evolution.
/// This function has no flag argument. Normalization shapes positions, never
/// speeds. `armc`/`lat_deg`/`eps_deg`
/// are degrees, `hsys` any selector form [`houses_armc`] accepts (the
/// folded character drives both the entry dispatch and the sign-pinned
/// speed rule), `ascmc9_deg` the Sun declination for `I`/`i` (ignored
/// otherwise). Returns 12 cusps (36 for Gauquelin `G`), eight ASCMC angles,
/// the same number of cusp speeds and eight angle speeds: Rust `Vec` + `[f64; 8]`
/// `HouseSpeedSolution`, positions `degnorm`-shaped, speeds in deg/day.
pub fn houses_armc_ex2(
    armc: f64,
    lat_deg: f64,
    eps_deg: f64,
    hsys: char,
    ascmc9_deg: f64,
) -> Result<HouseSpeedSolution> {
    let (cusps, ascmc) = houses_armc(armc, lat_deg, eps_deg, hsys, ascmc9_deg)?;
    let (cusps_before, ascmc_before) = houses_armc(
        armc - ARMC_SPEED_STEP_DEG,
        lat_deg,
        eps_deg,
        hsys,
        ascmc9_deg,
    )?;
    let (cusps_after, ascmc_after) = houses_armc(
        armc + ARMC_SPEED_STEP_DEG,
        lat_deg,
        eps_deg,
        hsys,
        ascmc9_deg,
    )?;
    let code = fold_hsys_case(hsys);
    let cusps_speed: Vec<f64> = cusps_after
        .iter()
        .zip(cusps_before.iter())
        .map(|(&after, &before)| armc_fd_speed_deg_per_day(after, before, true, code))
        .collect();
    let mut ascmc_speed = [0.0f64; N_ASCMC];
    for i in 0..N_ASCMC {
        ascmc_speed[i] = armc_fd_speed_deg_per_day(ascmc_after[i], ascmc_before[i], false, code);
    }
    Ok((cusps, ascmc, cusps_speed, ascmc_speed))
}

/// Plain centered-difference cusp rate in deg/day: `d/(2·_DT_DAYS)` with `_DT_DAYS = 1/4096` d ≈ 21 s
/// half-step — an exact binary fraction so `tjdut ± h` stays symmetric
/// at JD ~2.45e6, where a one-second step loses ~5 digits to roundoff
/// (ULP ≈ 46 µs), biasing every rate low by several arcsec/day.
///
/// Single-fold wrap (strict `> 180`/`< −180`); `> 90°` reports `0.0`;
/// sign-pinned precompute (once per call) with 30° ± 1e-6 reporting
/// `0.0` (the 30° step reads 30/(2·1/4096) = 61440 deg/day on all
/// twelve cusps). Arguments are degrees in (`after`/`before`,
/// precomputed sign-pinned flag), rate in deg/day out (UT1 days);
/// scalar on both sides.
pub fn ex2_cusp_rate_deg_per_day(after: f64, before: f64, sign_pinned: bool) -> f64 {
    let d = wrap_speed_diff_deg(after - before);
    if d.abs() > SPEED_DISCONTINUITY_DEG {
        return 0.0;
    }
    if sign_pinned && (d.abs() - SIGN_PINNED_STEP_DEG).abs() < SIGN_PINNED_STEP_TOL_DEG {
        return 0.0;
    }
    d / (2.0 * EX2_DT_DAYS)
}

/// Richardson-extrapolated angle rate in deg/day: wrapped central differences at coarse
/// `±h` and fine `±h/2 = ±1/8192` d — both exact binary fractions above
/// the roundoff floor — combined `D = [4·D(h/2) − D(h)]/3` to cancel the
/// leading `O(h²)` term (`O(h²) → O(h⁴)`).
///
/// Angles need this because they carry large `f'''` on the
/// latitude-amplified Ascendant/Vertex/co-ascendant slots: the fixed
/// 1/4096 truncation reaches a few arcsec/day at mid-latitude and
/// hundreds near the polar circle, while MC/ARMC sit at the roundoff
/// floor — and refining the step alone pushes the slow slots into
/// roundoff. Either-step `> 90°` reports `0.0` (no sign-pinned arm —
/// angles only). Arguments are degrees in (coarse/fine sample pairs),
/// rate in deg/day out; scalar on both sides.
pub fn ex2_angle_rate_richardson_deg_per_day(
    plus_h: f64,
    minus_h: f64,
    plus_q: f64,
    minus_q: f64,
) -> f64 {
    let d_h = wrap_speed_diff_deg(plus_h - minus_h);
    let d_q = wrap_speed_diff_deg(plus_q - minus_q);
    if d_h.abs() > SPEED_DISCONTINUITY_DEG || d_q.abs() > SPEED_DISCONTINUITY_DEG {
        return 0.0;
    }
    let d_coarse = d_h / (2.0 * EX2_DT_DAYS);
    let d_fine = d_q / (2.0 * EX2_DT_HALF_DAYS);
    (4.0 * d_fine - d_coarse) / 3.0
}

/// JD-based total-time-derivative speeds over an injected
/// degree-layer `houses_ex` sampler (provider-agnostic
/// contract; the JD acquisition stays with the epoch-entry layer): finite-differencing the
/// full `houses_ex` solution in time captures the ARMC rate, obliquity
/// drift, nutation and the `FLG_SIDEREAL` ayanamsha — the genuine
/// derivative of the reported positions for every system, including
/// iterative Placidus/Koch, which no analytic approximation reproduces.
///
/// `sample` maps a Julian Day (UT1) to the `(cusps, ascmc)` pair the
/// release `houses_ex(tjdut ± h, lat, lon, hsys, flags & ~FLG_RADIANS)`
/// calls produce — same latitude/longitude/system/degree flags at every
/// offset, five evaluations in the release order (center, `−h`, `+h`,
/// `−h/2`, `+h/2`; the fine pair feeds angles only). Cusps take the
/// plain centered rate, angles the Richardson rate; `folded_hsys`
/// drives the sign-pinned cusp arm (the caller folds the int code
/// at the public boundary; this kernel takes the folded character. Center
/// positions pass through unchanged. Speeds stay
/// deg/day (never radians); analytic Porphyry/Whole-Sign/Aries speeds
/// that are not the derivative of their own cusps are intentionally not
/// replicated. A refused sample fails the whole call — never a silent
/// zero. Returns `(cusps, auxiliary_angles, cusp_rates, angle_rates)` with
/// owned 12/36-element cusp vectors and eight-element angle arrays.
///
/// The callback must return the same cusp count and order at every sampled
/// epoch and supply degree-valued positions, with RADIANS removed by its caller.
/// This precondition is not checked: the two coarse cusp vectors are zipped,
/// so mismatched lengths truncate the rate vector even if the central vector
/// is longer. Fine samples contribute auxiliary angles only. Callback errors
/// abort without a partial result. Central positions and computed rates return
/// as owned values; no callback or additional sample history is retained.
pub fn houses_ex2_via<F>(sample: &F, tjdut: f64, folded_hsys: char) -> Result<HouseSpeedSolution>
where
    F: Fn(f64) -> Result<(Vec<f64>, [f64; N_ASCMC])>,
{
    houses_ex2_storage(sample, tjdut, folded_hsys)
}

/// Storage changes only ownership. Public Vec callbacks retain their arbitrary
/// lengths, including the original zip truncation, sample order and errors.
pub(super) fn houses_ex2_storage<F, C>(
    sample: &F,
    tjdut: f64,
    folded_hsys: char,
) -> Result<(C, [f64; N_ASCMC], C, [f64; N_ASCMC])>
where
    F: Fn(f64) -> Result<(C, [f64; N_ASCMC])>,
    C: AsRef<[f64]> + FromIterator<f64>,
{
    let (cusps, ascmc) = sample(tjdut)?;
    let (cusps_minus, ascmc_minus) = sample(tjdut - EX2_DT_DAYS)?;
    let (cusps_plus, ascmc_plus) = sample(tjdut + EX2_DT_DAYS)?;
    let (_, ascmc_hminus) = sample(tjdut - EX2_DT_HALF_DAYS)?;
    let (_, ascmc_hplus) = sample(tjdut + EX2_DT_HALF_DAYS)?;
    let sign_pinned = is_sign_pinned(folded_hsys);
    let cusps_speed: C = cusps_plus
        .as_ref()
        .iter()
        .zip(cusps_minus.as_ref().iter())
        .map(|(&after, &before)| ex2_cusp_rate_deg_per_day(after, before, sign_pinned))
        .collect();
    let mut ascmc_speed = [0.0f64; N_ASCMC];
    for i in 0..N_ASCMC {
        ascmc_speed[i] = ex2_angle_rate_richardson_deg_per_day(
            ascmc_plus[i],
            ascmc_minus[i],
            ascmc_hplus[i],
            ascmc_hminus[i],
        );
    }
    Ok((cusps, ascmc, cusps_speed, ascmc_speed))
}

#[cfg(test)]
mod exact_house_cache_tests {
    use super::*;

    fn words<C: IntoIterator<Item = f64>>((cusps, angles): (C, [f64; N_ASCMC])) -> Vec<u64> {
        cusps.into_iter().chain(angles).map(f64::to_bits).collect()
    }

    #[test]
    fn geometry_cache_preserves_system_site_obliquity_and_sun_inputs() {
        for system in "PKRCEAWOBTMXVHFUNSLQDIiJYGg?".chars() {
            for (armc, lat, eps, declination) in [
                (0.0, 41.9, 23.4, -23.5),
                (123.45, -30.0, 24.0, 10.0),
                (359.999, 0.0, 23.0, 0.0),
            ] {
                let expected =
                    words(houses_armc_uncached(armc, lat, eps, system, declination).unwrap());
                for _ in 0..3 {
                    assert_eq!(
                        words(
                            houses_armc_from_normalized(armc, lat, eps, system, declination)
                                .unwrap()
                        ),
                        expected
                    );
                }
            }
        }
        // A cached valid location cannot admit a later invalid or polar one.
        assert!(houses_armc(123.45, f64::NAN, 23.4, 'P', 0.0).is_err());
        assert!(matches!(
            houses_armc(123.45, 80.0, 23.4, 'P', 0.0),
            Err(Error::PolarCircle { .. })
        ));
    }
}
