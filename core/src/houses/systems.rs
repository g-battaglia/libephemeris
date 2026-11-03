//! Systems implementation for the parent module.
//! Public entry points and types are exposed through the parent API.

use super::*;

/// Native sine and cosine for an angle in radians.
#[inline]
pub(super) fn sin_cos(x: f64) -> (f64, f64) {
    x.sin_cos()
}

/// Normalize an angle to `[0, 360)`: worked inputs `45 / -45 / 370 / -370 /
/// 360 / 720` return `45 / 315 / 10 / 350 / 0 / 0`.
///
/// - Inputs with `|x| < 1e-13` (exclusive boundary) return exactly `0.0`
///   before the wrap, so a tiny negative never reaches the modulo. The
///   band sits above the modulo roundoff floor near zero (one ulp of 360
///   is ~5.7e-14) and far below angular resolution (`1e-13` deg =
///   `3.6e-10` arcsec).
/// - A modulo artifact of exactly `360.0` (roundoff can produce it for
///   tiny negatives such as `-1e-17`) snaps to `0.0`, keeping the
///   documented `[0, 360)` contract; legitimate near-360 values pass.
pub fn degnorm_deg(x: f64) -> f64 {
    crate::utils::degnorm(x)
}

/// Cardinal 13-slot cusp array (`_init_cardinal_cusps` ): index 0
/// unused, house 1 = Ascendant, house 10 = MC, houses 4/7 their
/// 180° opposites.
pub(super) fn init_cardinal_cusps_13(asc: f64, mc: f64) -> [f64; 13] {
    let mut cusps = [0.0_f64; 13];
    cusps[1] = asc;
    cusps[10] = mc;
    cusps[7] = wrap360(asc + 180.0);
    cusps[4] = wrap360(mc + 180.0);
    cusps
}

/// Opposite-house assignment (`_set_opposite_cusps` ): houses
/// 5/6/8/9 are the 180° opposites of 11/12/2/3.
pub(super) fn set_opposite_cusps_13(cusps: &mut [f64; 13]) {
    cusps[5] = wrap360(cusps[11] + 180.0);
    cusps[6] = wrap360(cusps[12] + 180.0);
    cusps[8] = wrap360(cusps[2] + 180.0);
    cusps[9] = wrap360(cusps[3] + 180.0);
}

/// Strip the unused index-0 slot to the public 12-cusp order (houses
/// 1–12), mirroring the `cusps[1:13]` slice of the release entry.
/// Per-cusp `degnorm` shaping is the separate [`degnorm_deg`] step at
/// the public entry, as in the release `houses` return.
pub(super) fn strip_index0(cusps: [f64; 13]) -> [f64; 12] {
    let mut out = [0.0_f64; 12];
    out.copy_from_slice(&cusps[1..13]);
    out
}

/// Equal houses: the simplest
/// system, 30° ecliptic divisions from house 1 at the Ascendant.
///
/// `λᵢ = (λ_Asc + (i−1)·30°) mod 360°` for `i = 1..12`. Works at all
/// latitudes including the polar circle; the MC need not coincide with
/// house 10; computationally exact (no astronomical calculation beyond
/// the wrap). Returns `[f64; 12]` houses 1–12.
///
/// `asc`: ecliptic longitude of the Ascendant in degrees, finite;
/// non-finite inputs propagate release-identically (see the module
/// contract above).
pub fn houses_equal(asc: f64) -> [f64; 12] {
    let mut cusps = [0.0_f64; 13];
    for (i, slot) in cusps.iter_mut().enumerate().skip(1) {
        *slot = wrap360(asc + (i as f64 - 1.0) * 30.0);
    }
    strip_index0(cusps)
}

/// Whole-sign houses: the
/// Hellenistic system, one zodiac sign per house with house 1 at 0° of
/// the Ascendant's sign (Holden 1977).
///
/// `start = floor((λ_Asc + 1e-9°) / 30°)·30°`, then
/// `λᵢ = (start + (i−1)·30°) mod 360°` for `i = 1..12`. The `1e-9°`
/// nudge ([`WHOLE_SIGN_NUDGE_DEG`]) resolves an Ascendant within one float
/// ULP below an exact sign cusp (e.g. `180 − 1.4e-13` at cardinal ARMC)
/// to the new sign, as the shipped code does. E.g. Ascendant 15° Taurus
/// (45°) starts house 1 at 30°. Works at all latitudes; cusps always sit
/// at sign 0°. Returns `[f64; 12]` houses 1–12.
///
/// `asc`: ecliptic longitude of the Ascendant in degrees, finite;
/// non-finite inputs propagate release-identically (see the module
/// contract above).
pub fn houses_whole_sign(asc: f64) -> [f64; 12] {
    let mut cusps = [0.0_f64; 13];
    let start = ((asc + WHOLE_SIGN_NUDGE_DEG) / 30.0).floor() * 30.0;
    for (i, slot) in cusps.iter_mut().enumerate().skip(1) {
        *slot = wrap360(start + (i as f64 - 1.0) * 30.0);
    }
    strip_index0(cusps)
}

/// Vehlow equal houses: the Equal
/// variant with the Ascendant 15° into house 1 (Johannes Vehlow,
/// *Erlebte Sternenwelt*).
///
/// `start = (λ_Asc − 15°) mod 360°`, then
/// `λᵢ = (start + (i−1)·30°) mod 360°` for `i = 1..12`. The Ascendant is
/// exactly 15° into house 1; works at all latitudes including polar;
/// computationally exact. Returns `[f64; 12]` houses 1–12.
///
/// `asc`: ecliptic longitude of the Ascendant in degrees, finite;
/// non-finite inputs propagate release-identically (see the module
/// contract above).
pub fn houses_vehlow(asc: f64) -> [f64; 12] {
    let mut cusps = [0.0_f64; 13];
    let start = wrap360(asc - VEHLOW_OFFSET_DEG);
    for (i, slot) in cusps.iter_mut().enumerate().skip(1) {
        *slot = wrap360(start + (i as f64 - 1.0) * 30.0);
    }
    strip_index0(cusps)
}

/// Equal houses from MC: 30° houses
/// with house 10 at the MC (`H10 = MC`, `H1 = MC − 270° = MC + 90°`,
/// `λᵢ = (MC + (i−10)·30°) mod 360°`). The Ascendant is not a cusp in
/// this system (kept as an unused parameter for the release call shape;
/// the Ascendant lives in the ASCMC vector). Works at all latitudes.
/// Returns `[f64; 12]` houses 1–12.
///
/// `mc`: ecliptic longitude of the MC in degrees, finite (`asc` is
/// unused, kept for the release call shape); non-finite inputs
/// propagate release-identically (see the module contract above).
pub fn houses_equal_mc(_asc: f64, mc: f64) -> [f64; 12] {
    let mut cusps = [0.0_f64; 13];
    for (i, slot) in cusps.iter_mut().enumerate().skip(1) {
        *slot = wrap360(mc + (i as f64 - 10.0) * 30.0);
    }
    strip_index0(cusps)
}

/// Porphyry houses: space-based
/// trisection of each Ascendant/MC quadrant along the ecliptic — simple,
/// no trigonometry — and the default fallback system
/// ([`DEFAULT_FALLBACK_HSYS`]) when Placidus/Koch fail.
///
/// With `λ_Asc` the Ascendant and `λ_MC` the MC: the MC→Asc step is
/// `((λ_Asc − λ_MC) mod 360°)/3`, giving `λ₁₁ = λ_MC + step` and
/// `λ₁₂ = λ_MC + 2·step`; the Asc→IC step is
/// `((λ_IC − λ_Asc) mod 360°)/3` with `λ_IC = λ_MC + 180°`, giving
/// `λ₂ = λ_Asc + step` and `λ₃ = λ_Asc + 2·step`; houses 5/6/8/9 are the
/// 180° opposites of 11/12/2/3. Works at all latitudes including polar.
/// Returns `[f64; 12]` houses 1–12.
///
/// `asc`/`mc`: ecliptic longitudes of the Ascendant and MC in degrees,
/// finite; non-finite inputs propagate release-identically (see the
/// module contract above).
pub fn houses_porphyry(asc: f64, mc: f64) -> [f64; 12] {
    let mut cusps = init_cardinal_cusps_13(asc, mc);
    // Arc MC → Asc (houses 11, 12).
    let diff = wrap360(asc - mc);
    let step = diff / 3.0;
    cusps[11] = wrap360(mc + step);
    cusps[12] = wrap360(mc + 2.0 * step);
    // Arc Asc → IC (houses 2, 3).
    let ic = cusps[4];
    let diff = wrap360(ic - asc);
    let step = diff / 3.0;
    cusps[2] = wrap360(asc + step);
    cusps[3] = wrap360(asc + 2.0 * step);
    set_opposite_cusps_13(&mut cusps);
    strip_index0(cusps)
}

/// Sripati houses: midpoints of the
/// Porphyry cusps (Jyotiṣa bhāva-cakra, Śrīpati; Holden 1977) — each
/// Porphyry cusp becomes the center of the corresponding Sripati house.
///
/// `λ'ᵢ = λᵢ₋₁ + ((λᵢ − λᵢ₋₁) mod 360°)/2`, with house 1 midpointing
/// houses 12 → 1. Inherits Porphyry polar stability; Jyotish use; simple.
/// Returns `[f64; 12]` houses 1–12.
///
/// `asc`/`mc`: ecliptic longitudes of the Ascendant and MC in degrees,
/// finite; non-finite inputs propagate release-identically (see the
/// module contract above).
pub fn houses_sripati(asc: f64, mc: f64) -> [f64; 12] {
    let porphyry = houses_porphyry(asc, mc);
    // Re-expand to the 13-slot layout the release midpoint loop reads.
    let mut base = [0.0_f64; 13];
    base[1..13].copy_from_slice(&porphyry);
    let mut cusps = [0.0_f64; 13];
    for i in 1..13 {
        let prev = if i == 1 { 12 } else { i - 1 };
        let diff = wrap360(base[i] - base[prev]);
        cusps[i] = wrap360(base[prev] + diff / 2.0);
    }
    strip_index0(cusps)
}

/// Pullen SD house widths for one quadrant (`calc_house_sizes` ):
/// `(side, middle, side)` with `middle = 30 + d/2`, `side = 30 + d/4`
/// for the deviation `d = quadrant − 90`; a negative middle clamps to
/// zero and splits the quadrant evenly across the side houses.
///
/// `quadrant_deg`: quadrant span in degrees, finite; non-finite inputs
/// propagate release-identically (see the module contract above).
pub(super) fn pullen_sd_quadrant_sizes(quadrant_deg: f64) -> (f64, f64, f64) {
    let d = quadrant_deg - 90.0;
    let middle = 30.0 + d / 2.0;
    if middle < 0.0 {
        (quadrant_deg / 2.0, 0.0, quadrant_deg / 2.0)
    } else {
        let side = 30.0 + d / 4.0;
        (side, middle, side)
    }
}

/// Pullen SD houses: sinusoidal-delta
/// (neo-Porphyry) quadrant weighting in the release quadrant layout,
/// accumulated from each quadrant root (Walter Pullen, 1994; Astrolog docs).
///
/// For a quadrant of `q` degrees with deviation `d = q − 90°`, the three
/// house widths follow the pattern `(x + n, x + 2n, x + n)` with `x = 30`
/// and `n = d/4`: sides `30 + d/4`, middle `30 + d/2` (a negative middle
/// clamps to `0`, splitting the remainder as `q/2` per side). Ecliptic-only;
/// works at all latitudes; larger quadrants get larger houses with the
/// middle house absorbing the variation. Quadrant→house map: Q1 (from MC)
/// fills 11/12, Q2 (from Ascendant) 2/3, Q3 (from IC) 5/6, Q4 (from
/// Descendant) 8/9. Returns `[f64; 12]` houses 1–12.
///
/// `asc`/`mc`: ecliptic longitudes of the Ascendant and MC in degrees,
/// finite; non-finite inputs propagate release-identically (see the
/// module contract above).
pub fn houses_pullen_sd(asc: f64, mc: f64) -> [f64; 12] {
    let mut cusps = init_cardinal_cusps_13(asc, mc);
    let ic = cusps[4];
    let desc = cusps[7];
    let q1 = wrap360(asc - mc);
    let q2 = wrap360(ic - asc);
    let q3 = wrap360(desc - ic);
    let q4 = wrap360(mc - desc);
    let s1 = pullen_sd_quadrant_sizes(q1);
    cusps[11] = wrap360(mc + s1.0);
    cusps[12] = wrap360(cusps[11] + s1.1);
    let s2 = pullen_sd_quadrant_sizes(q2);
    cusps[2] = wrap360(asc + s2.0);
    cusps[3] = wrap360(cusps[2] + s2.1);
    let s3 = pullen_sd_quadrant_sizes(q3);
    cusps[5] = wrap360(ic + s3.0);
    cusps[6] = wrap360(cusps[5] + s3.1);
    let s4 = pullen_sd_quadrant_sizes(q4);
    cusps[8] = wrap360(desc + s4.0);
    cusps[9] = wrap360(cusps[8] + s4.1);
    strip_index0(cusps)
}

/// Pullen SR house widths for one quadrant (same operation order as
/// `calc_house_sizes_ratio`): `(side, middle, side) = (r·x, x, r·x)` with the ratio `r`
/// solving `r³(2 + r) / (2r + 1) = (180 − q) / q` by Newton–Raphson from
/// `r = 1` (same iteration budget, derivative floor, positivity
/// clamp and stop rule); a quadrant within
/// [`PULLEN_SR_EQUAL_BAND_DEG`] of 90° takes three exact 30° houses.
///
/// `quadrant_deg`: quadrant span in degrees, finite; non-finite inputs
/// propagate release-identically (see the module contract above).
pub(super) fn pullen_sr_quadrant_sizes(quadrant_deg: f64) -> (f64, f64, f64) {
    let q = quadrant_deg;
    let q_opp = 180.0 - q;
    if (q - 90.0).abs() < PULLEN_SR_EQUAL_BAND_DEG {
        return (30.0, 30.0, 30.0);
    }
    let target = if q > PULLEN_SR_DIVISOR_GUARD {
        q_opp / q
    } else {
        PULLEN_SR_EMPTY_QUADRANT_TARGET
    };
    let mut r = 1.0_f64;
    for _ in 0..PULLEN_SR_MAX_ITER {
        let numerator = r * r * r * (2.0 + r);
        let denominator = 2.0 * r + 1.0;
        let f = numerator / denominator - target;
        let dn = 3.0 * r * r * (2.0 + r) + r * r * r;
        let dd = 2.0;
        let df = (dn * denominator - numerator * dd) / (denominator * denominator);
        if df.abs() < PULLEN_SR_DERIVATIVE_FLOOR {
            break;
        }
        let mut r_new = r - f / df;
        if r_new < PULLEN_SR_R_MIN {
            r_new = PULLEN_SR_R_MIN;
        }
        if (r_new - r).abs() < PULLEN_SR_CONVERGENCE {
            break;
        }
        r = r_new;
    }
    let x = if 2.0 * r + 1.0 > PULLEN_SR_DIVISOR_GUARD {
        q / (2.0 * r + 1.0)
    } else {
        q / 3.0
    };
    (r * x, x, r * x)
}

/// Pullen SR houses: sinusoidal-ratio
/// quadrant weighting (2016 improvement over SD; Astrolog docs), accumulated
/// in the same quadrant layout as [`houses_pullen_sd`] but with
/// `pullen_sr_quadrant_sizes` widths.
///
/// For a quadrant of `q` degrees against the opposite `180 − q`, the small
/// side solves `(rx, x, rx)` summing to `q` and the large side
/// `(r³x, r⁴x, r³x)` summing to `180 − q`, joined by the dividing equation
/// `r³(2 + r)/(2r + 1) = (180 − q)/q` (solved numerically); ratio
/// multipliers keep every size non-negative, unlike the SD form. Works at
/// all latitudes. Returns `[f64; 12]` houses 1–12.
///
/// `asc`/`mc`: ecliptic longitudes of the Ascendant and MC in degrees,
/// finite; non-finite inputs propagate release-identically (see the
/// module contract above).
pub fn houses_pullen_sr(asc: f64, mc: f64) -> [f64; 12] {
    let mut cusps = init_cardinal_cusps_13(asc, mc);
    let ic = cusps[4];
    let desc = cusps[7];
    let q1 = wrap360(asc - mc);
    let q2 = wrap360(ic - asc);
    let q3 = wrap360(desc - ic);
    let q4 = wrap360(mc - desc);
    let s1 = pullen_sr_quadrant_sizes(q1);
    cusps[11] = wrap360(mc + s1.0);
    cusps[12] = wrap360(cusps[11] + s1.1);
    let s2 = pullen_sr_quadrant_sizes(q2);
    cusps[2] = wrap360(asc + s2.0);
    cusps[3] = wrap360(cusps[2] + s2.1);
    let s3 = pullen_sr_quadrant_sizes(q3);
    cusps[5] = wrap360(ic + s3.0);
    cusps[6] = wrap360(cusps[5] + s3.1);
    let s4 = pullen_sr_quadrant_sizes(q4);
    cusps[8] = wrap360(desc + s4.0);
    cusps[9] = wrap360(cusps[8] + s4.1);
    strip_index0(cusps)
}

/// Bound a direction cosine to [-1,1], preserving NaN and signed zero.
#[inline]
pub(super) fn clamp_unit(x: f64) -> f64 {
    x.clamp(-1.0, 1.0)
}

/// Right ascension to ecliptic longitude (Smart/Meeus spherical-trigonometry
/// construction): `tan(λ) = sin(α) / (cos(ε)·cos(α) − sin(ε)·tan(φ))`
/// with `α` the right ascension, `ε` the obliquity, `φ` the pole height,
/// and `atan2` quadrant resolution, covering the polar pole-height
/// degeneracies, the 0/0 continuous limit and the cardinal snaps.
///
/// `ra_deg` is degrees; `pole_height_deg` is degrees; the obliquity
/// travels pre-computed (`sin_obliquity`, `cos_obliquity`) so the
/// Campanus caller evaluates `sin`/`cos` once for its four cusps.
/// Returns the longitude in `[0, 360)`; scalar on both sides (no cusp
/// layout applies here).
pub fn ra_to_ecliptic_longitude(
    ra_deg: f64,
    pole_height_deg: f64,
    sin_obliquity: f64,
    cos_obliquity: f64,
) -> f64 {
    if (90.0 - pole_height_deg).abs() < RA_LON_NEAR_ZERO_DEG {
        return 180.0;
    }
    if (90.0 + pole_height_deg).abs() < RA_LON_NEAR_ZERO_DEG {
        return 0.0;
    }
    let ra_rad = wrap360(ra_deg).to_radians();
    let tan_pole = pole_height_deg.to_radians().tan();
    let (sin_ra, cos_ra) = sin_cos(ra_rad);
    let numerator = sin_ra;
    let denominator = cos_obliquity * cos_ra - sin_obliquity * tan_pole;
    if numerator.abs() < RA_LON_SINGULARITY_FLOOR && denominator.abs() < RA_LON_SINGULARITY_FLOOR {
        return wrap360(ra_deg);
    }
    let longitude = wrap360(numerator.atan2(denominator).to_degrees());
    for cardinal in [90.0, 180.0, 270.0] {
        if (longitude - cardinal).abs() < RA_LON_NEAR_ZERO_DEG {
            return cardinal;
        }
    }
    if (longitude - 360.0).abs() < RA_LON_NEAR_ZERO_DEG || longitude.abs() < RA_LON_NEAR_ZERO_DEG {
        return 0.0;
    }
    longitude
}

/// Regiomontanus houses:
/// equatorial 30° arcs projected through the celestial poles
/// (Regiomontanus, *Tabulae directionum*, 1467, printed 1490; Holden
/// 1977).
///
/// With `H = 30°/60°/120°/150°` selecting houses 11/12/2/3:
/// `tan(P) = tan(φ)·sin(H)` with `R = ARMC + H − 90°` and
/// `λ = atan2(cos R, −(sin R·cos ε + tan P·sin ε))`; houses 5/6/8/9 are
/// the 180° opposites. Arguments are degrees: the active `armc_deg` (the
/// entry prologue owns the raising-selector polar-circle refusal and the
/// R/C/T/J MC-flip; this kernel takes the active ARMC plus the
/// precomputed Ascendant/MC cardinals as given), latitude, true
/// obliquity, and the `asc_deg`/`mc_deg` cardinals seeding houses 1/10.
/// Returns `[f64; 12]` houses 1–12.
///
/// Per-cusp `degnorm` shaping is the separate [`degnorm_deg`] step at
/// the public entry (see the module contract above).
pub fn houses_regiomontanus(
    armc_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
    asc_deg: f64,
    mc_deg: f64,
) -> [f64; 12] {
    fn calc_cusp(armc_deg: f64, rad_lat: f64, rad_eps: f64, offset_deg: f64) -> f64 {
        let h_rad = offset_deg.to_radians();
        let tan_pole = rad_lat.tan() * h_rad.sin();
        let r_deg = wrap360(armc_deg + offset_deg - 90.0);
        let r_rad = r_deg.to_radians();
        let (sin_r, cos_r) = sin_cos(r_rad);
        let (sin_eps, cos_eps) = sin_cos(rad_eps);
        let num = cos_r;
        let den = -(sin_r * cos_eps + tan_pole * sin_eps);
        wrap360(num.atan2(den).to_degrees())
    }

    let mut cusps = init_cardinal_cusps_13(asc_deg, mc_deg);
    let rad_lat = lat_deg.to_radians();
    let rad_eps = eps_deg.to_radians();
    cusps[11] = calc_cusp(armc_deg, rad_lat, rad_eps, 30.0);
    cusps[12] = calc_cusp(armc_deg, rad_lat, rad_eps, 60.0);
    cusps[2] = calc_cusp(armc_deg, rad_lat, rad_eps, 120.0);
    cusps[3] = calc_cusp(armc_deg, rad_lat, rad_eps, 150.0);
    set_opposite_cusps_13(&mut cusps);
    strip_index0(cusps)
}

/// Campanus houses: prime-vertical
/// 30° arcs projected through the horizon north/south points (Campanus
/// of Novara, 13th century; Holden 1977; pole geometry after Smart
/// Ch. 3 / Meeus Ch. 13).
///
/// The north point is the pole of the prime vertical: each house circle
/// through the north/Desc/south points, perpendicular to the prime
/// vertical at the PV offset `h`, has pole height
/// `arcsin(sin φ·cos h)` and pole hour angle
/// `atan2(sin h, cos h·cos φ)`; the equator crossing at
/// `ARMC + 90° − H_pole` converts via
/// [`ra_to_ecliptic_longitude`]. Cusps 12/11 sit at `+30°/+60°` (above
/// the horizon, between the Ascendant and the MC), cusps 2/3 at
/// `−30°/−60°` (below the horizon, between the Ascendant and the IC);
/// houses 5/6/8/9 are the 180° opposites. Arguments are degrees
/// (active ARMC, latitude, true obliquity, precomputed Ascendant/MC).
/// Returns `[f64; 12]` houses 1–12.
///
/// Per-cusp `degnorm` shaping is the separate [`degnorm_deg`] step at
/// the public entry (see the module contract above).
pub fn houses_campanus(
    armc_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
    asc_deg: f64,
    mc_deg: f64,
) -> [f64; 12] {
    fn cusp_from_pv_offset(
        armc_deg: f64,
        sin_lat: f64,
        cos_lat: f64,
        sin_eps: f64,
        cos_eps: f64,
        h_deg: f64,
    ) -> f64 {
        let h = h_deg.to_radians();
        let (sin_h, cos_h) = sin_cos(h);
        let pole_height = clamp_unit(sin_lat * cos_h).asin().to_degrees();
        let h_pole = sin_h.atan2(cos_h * cos_lat).to_degrees();
        let ra_crossing = armc_deg + 90.0 - h_pole;
        ra_to_ecliptic_longitude(ra_crossing, pole_height, sin_eps, cos_eps)
    }

    let mut cusps = init_cardinal_cusps_13(asc_deg, mc_deg);
    let lat_rad = lat_deg.to_radians();
    let (sin_lat, cos_lat) = sin_cos(lat_rad);
    let eps_rad = eps_deg.to_radians();
    let (sin_eps, cos_eps) = sin_cos(eps_rad);
    cusps[12] = cusp_from_pv_offset(armc_deg, sin_lat, cos_lat, sin_eps, cos_eps, 30.0);
    cusps[11] = cusp_from_pv_offset(armc_deg, sin_lat, cos_lat, sin_eps, cos_eps, 60.0);
    cusps[2] = cusp_from_pv_offset(armc_deg, sin_lat, cos_lat, sin_eps, cos_eps, -30.0);
    cusps[3] = cusp_from_pv_offset(armc_deg, sin_lat, cos_lat, sin_eps, cos_eps, -60.0);
    set_opposite_cusps_13(&mut cusps);
    strip_index0(cusps)
}

/// Alcabitius houses: RA-arc
/// trisection MC→Ascendant (houses 11/12) and Ascendant→IC (houses 2/3)
/// with hour-circle projection (al-Qabīṣī, 10th century;
/// Burnett/Yamamoto/Yano 2004).
///
/// Only the Ascendant's right ascension is computed,
/// `RA_Asc = atan2(sin λ_Asc·cos ε, cos λ_Asc)`; the MC side uses the
/// active ARMC directly as the equatorial origin:
/// `arc = (RA_Asc − ARMC) mod 360°`, `step = arc/3`,
/// `RA₁₁/₁₂ = ARMC + step/2·step`, and with `RA_IC = ARMC + 180°`,
/// `arc2`/`step2` give `RA₂/₃ = RA_Asc + step2/2·step2`; each projects
/// back with `λ = atan2(sin RA, cos RA·cos ε)` (so house 11 sits `SA/3`
/// past the RAMC with `SA = 90° + AD_asc`). Houses 5/6/8/9 are the 180°
/// opposites. The release evaluates `math.radians(lat)` and discards
/// it; the latitude plays no further role in this construction.
/// Arguments are degrees (active ARMC, latitude unused, true obliquity,
/// precomputed Ascendant/MC). Returns `[f64; 12]` houses 1–12 (the
/// release 13-list with index 0 unused, rendered in Rust order).
///
/// Per-cusp `degnorm` shaping is the separate [`degnorm_deg`] step at
/// the public entry (see the module contract above).
pub fn houses_alcabitius(
    armc_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
    asc_deg: f64,
    mc_deg: f64,
) -> [f64; 12] {
    fn lon_from_ra(rad_eps: f64, ra_deg: f64) -> f64 {
        let ra_rad = ra_deg.to_radians();
        let (sin_ra, cos_ra) = sin_cos(ra_rad);
        let y = sin_ra;
        let x = cos_ra * rad_eps.cos();
        wrap360(y.atan2(x).to_degrees())
    }

    let mut cusps = init_cardinal_cusps_13(asc_deg, mc_deg);
    let rad_eps = eps_deg.to_radians();
    let _ = lat_deg;
    let asc_rad = asc_deg.to_radians();
    let (sin_asc, cos_asc) = sin_cos(asc_rad);
    let ra_asc = wrap360((sin_asc * rad_eps.cos()).atan2(cos_asc).to_degrees());
    let arc = wrap360(ra_asc - armc_deg);
    let step = arc / 3.0;
    cusps[11] = lon_from_ra(rad_eps, armc_deg + step);
    cusps[12] = lon_from_ra(rad_eps, armc_deg + 2.0 * step);
    let ra_ic = wrap360(armc_deg + 180.0);
    let arc2 = wrap360(ra_ic - ra_asc);
    let step2 = arc2 / 3.0;
    cusps[2] = lon_from_ra(rad_eps, ra_asc + step2);
    cusps[3] = lon_from_ra(rad_eps, ra_asc + 2.0 * step2);
    set_opposite_cusps_13(&mut cusps);
    strip_index0(cusps)
}

/// Morinus houses: equatorial 30°
/// sections projected through the ecliptic poles (Morin, *Astrologia
/// Gallica*, 1661, Book XVII).
///
/// `RA = ARMC + (i−10)·30°`, `tan λ = tan RA·cos ε`,
/// `λ = atan2(sin RA·cos ε, cos RA)`. Location-independent: the
/// latitude, Ascendant and MC only shape the call, never the cusps —
/// house 1 is the `ARMC + 90°` projection, not the horizon Ascendant,
/// and house 10 is the `ARMC` projection, not the MC. Houses
/// 10/11/12/1/2/3 fill directly with 4/5/6/7/8/9 as their 180°
/// opposites (own opposites, not the shared helper). Arguments are
/// degrees (active ARMC, latitude/Ascendant/MC unused, true obliquity).
/// Returns `[f64; 12]` houses 1–12.
///
/// Per-cusp `degnorm` shaping is the separate [`degnorm_deg`] step at
/// the public entry (see the module contract above).
pub fn houses_morinus(
    armc_deg: f64,
    _lat_deg: f64,
    eps_deg: f64,
    _asc_deg: f64,
    _mc_deg: f64,
) -> [f64; 12] {
    fn lon_from_ra(rad_eps: f64, ra_deg: f64) -> f64 {
        let ra_rad = ra_deg.to_radians();
        let (sin_ra, cos_ra) = sin_cos(ra_rad);
        let y = sin_ra * rad_eps.cos();
        let x = cos_ra;
        wrap360(y.atan2(x).to_degrees())
    }

    let rad_eps = eps_deg.to_radians();
    let mut cusps = [0.0_f64; 13];
    cusps[10] = lon_from_ra(rad_eps, armc_deg);
    cusps[11] = lon_from_ra(rad_eps, armc_deg + 30.0);
    cusps[12] = lon_from_ra(rad_eps, armc_deg + 60.0);
    cusps[1] = lon_from_ra(rad_eps, armc_deg + 90.0);
    cusps[2] = lon_from_ra(rad_eps, armc_deg + 120.0);
    cusps[3] = lon_from_ra(rad_eps, armc_deg + 150.0);
    cusps[4] = wrap360(cusps[10] + 180.0);
    cusps[5] = wrap360(cusps[11] + 180.0);
    cusps[6] = wrap360(cusps[12] + 180.0);
    cusps[7] = wrap360(cusps[1] + 180.0);
    cusps[8] = wrap360(cusps[2] + 180.0);
    cusps[9] = wrap360(cusps[3] + 180.0);
    strip_index0(cusps)
}

/// Meridian houses: meridional
/// 30° right-ascension intervals from the MC projected through the
/// celestial poles (Zariel/David Cope, c. 1910; Holden 1977).
///
/// `RA = ARMC + (i−10)·30°`, `tan λ = tan RA / cos ε`,
/// `λ = atan2(sin RA, cos RA·cos ε)` — contrast Morinus
/// `atan2(sin RA·cos ε, cos RA)` with Meridian
/// `atan2(sin RA, cos RA·cos ε)`: the same direct-fill layout as
/// [`houses_morinus`] (10/11/12/1/2/3 direct, 4–9 as 180° opposites)
/// with the transposed projection. Location- and angle-independent:
/// latitude, Ascendant and MC only shape the call (degrees, with the
/// active ARMC and true obliquity). Returns `[f64; 12]` houses 1–12
///
/// Per-cusp `degnorm` shaping is the separate [`degnorm_deg`] step at
/// the public entry (see the module contract above).
pub fn houses_meridian(
    armc_deg: f64,
    _lat_deg: f64,
    eps_deg: f64,
    _asc_deg: f64,
    _mc_deg: f64,
) -> [f64; 12] {
    fn lon_from_ra(rad_eps: f64, ra_deg: f64) -> f64 {
        let ra_rad = ra_deg.to_radians();
        let (sin_ra, cos_ra) = sin_cos(ra_rad);
        let y = sin_ra;
        let x = cos_ra * rad_eps.cos();
        wrap360(y.atan2(x).to_degrees())
    }

    let rad_eps = eps_deg.to_radians();
    let mut cusps = [0.0_f64; 13];
    cusps[10] = lon_from_ra(rad_eps, armc_deg);
    cusps[11] = lon_from_ra(rad_eps, armc_deg + 30.0);
    cusps[12] = lon_from_ra(rad_eps, armc_deg + 60.0);
    cusps[1] = lon_from_ra(rad_eps, armc_deg + 90.0);
    cusps[2] = lon_from_ra(rad_eps, armc_deg + 120.0);
    cusps[3] = lon_from_ra(rad_eps, armc_deg + 150.0);
    cusps[4] = wrap360(cusps[10] + 180.0);
    cusps[5] = wrap360(cusps[11] + 180.0);
    cusps[6] = wrap360(cusps[12] + 180.0);
    cusps[7] = wrap360(cusps[1] + 180.0);
    cusps[8] = wrap360(cusps[2] + 180.0);
    cusps[9] = wrap360(cusps[3] + 180.0);
    strip_index0(cusps)
}

/// Bound the ascensional product `tan(latitude)*tan(declination)` to
/// [-1,1] before inverse sine. Non-finite NaN geometry remains unknown.
#[inline]
pub(super) fn clamp_ascensional_product(prod: f64) -> f64 {
    prod.clamp(-1.0, 1.0)
}

/// Solve one Placidus cusp by fixed-point iteration.
///
/// `offset_deg` selects the cusp and its initial right-ascension displacement:
/// 30/60 select houses 11/12 from the upper meridian, and 120/150 select
/// houses 2/3 from the lower meridian. Latitude and obliquity are radians;
/// ARMC, offsets and the returned ecliptic longitude are degrees.
///
/// Each pass recomputes `tan δ = sin RA·tan ε`, the clamped ascensional
/// difference, and the next `RA`, stopping at `|ΔRA| < 1e-7°`
/// ([`PLACIDUS_CONVERGENCE_DEG`]); typical convergence takes ~15–20
/// iterations. The clamp keeps iterating through near-polar candidates
/// whose mid-iteration `|tan φ·tan δ|` exceeds 1 even though the
/// converged cusp is well-defined. Returns `None` when the iteration
/// diverges or exhausts
/// [`PLACIDUS_MAX_ITERATIONS`], letting the caller fall back to
/// [`bisect_placidus`].
fn iterate_placidus(armc_deg: f64, rad_lat: f64, rad_eps: f64, offset_deg: f64) -> Option<f64> {
    let mut ra = wrap360(armc_deg + offset_deg);
    let mut prev_diff = f64::INFINITY;
    for iteration in 0..PLACIDUS_MAX_ITERATIONS {
        let sin_ra = ra.to_radians().sin();
        let tan_dec = sin_ra * rad_eps.tan();
        let prod = clamp_ascensional_product(rad_lat.tan() * tan_dec);
        let ad_deg = prod.asin().to_degrees();
        let mut new_ra = ra;
        if offset_deg == 30.0 {
            let h_deg = (90.0 + ad_deg) / 3.0;
            new_ra = wrap360(armc_deg + h_deg);
        } else if offset_deg == 60.0 {
            let h_deg = 2.0 * (90.0 + ad_deg) / 3.0;
            new_ra = wrap360(armc_deg + h_deg);
        } else if offset_deg == 120.0 {
            let h_deg = 2.0 * (90.0 - ad_deg) / 3.0;
            new_ra = wrap360(armc_deg + 180.0 - h_deg);
        } else if offset_deg == 150.0 {
            // House 3 uses one third of the nocturnal semi-arc.
            let h_deg = 1.0 * (90.0 - ad_deg) / 3.0;
            new_ra = wrap360(armc_deg + 180.0 - h_deg);
        }
        let mut diff = (new_ra - ra).abs();
        if diff > 180.0 {
            diff = 360.0 - diff;
        }
        ra = new_ra;
        if diff < PLACIDUS_CONVERGENCE_DEG {
            let ra_rad = ra.to_radians();
            let (sin_ra, cos_ra) = sin_cos(ra_rad);
            let y = sin_ra;
            let x = cos_ra * rad_eps.cos();
            return Some(wrap360(y.atan2(x).to_degrees()));
        }
        if iteration > PLACIDUS_DIVERGENCE_WARMUP_ITERATIONS
            && diff > prev_diff * PLACIDUS_DIVERGENCE_FACTOR
        {
            return None;
        }
        prev_diff = diff;
    }
    None
}

/// Ascensional difference in degrees at a degree-valued right ascension.
/// Latitude and obliquity are radians. The circumpolar product is clamped
/// consistently with the fixed-point solver, giving a result in `[-90, 90]`.
fn placidus_ad_deg_at(armc_term: f64, rad_lat: f64, rad_eps: f64) -> f64 {
    let tan_dec = armc_term.to_radians().sin() * rad_eps.tan();
    let prod = clamp_ascensional_product(rad_lat.tan() * tan_dec);
    prod.asin().to_degrees()
}

/// Placidus bisection fallback: the four cusp equations as continuous roots —
/// H11: `x = (90 + AD(ARMC + x))/3` on `[0, 60]`;
/// H12: `x = 2(90 + AD)/3` on `[0, 120]`;
/// H2: `x = 2(90 − AD(ARMC + 180 − x))/3` on `[0, 120]`;
/// H3: `x = (90 − AD)/3` on `[0, 60]` — i.e. `g(x) = x − frac·(90 ± AD)`
/// with these brackets halved [`PLACIDUS_BISECTION_ITERATIONS`] times.
///
/// Bisection exists because the fixed point oscillates within float
/// noise of the polar circle (e.g. ±66.55° at ε 23.4393°); the continuous
/// equation retains a root with a
/// sign change. Returns `None` only when the bracket itself fails
/// (unreachable with the clamped AD).
/// ARMC, the cusp selector and the returned longitude are degrees;
/// latitude and obliquity are radians, as in [`iterate_placidus`].
fn bisect_placidus(armc_deg: f64, rad_lat: f64, rad_eps: f64, offset_deg: f64) -> Option<f64> {
    let (frac, hi, from_ic) = if offset_deg == 30.0 {
        (1.0 / 3.0, 60.0, false)
    } else if offset_deg == 60.0 {
        (2.0 / 3.0, 120.0, false)
    } else if offset_deg == 120.0 {
        (2.0 / 3.0, 120.0, true)
    } else {
        (1.0 / 3.0, 60.0, true)
    };
    let g = |x: f64| -> f64 {
        if from_ic {
            let ad = placidus_ad_deg_at(armc_deg + 180.0 - x, rad_lat, rad_eps);
            x - frac * (90.0 - ad)
        } else {
            let ad = placidus_ad_deg_at(armc_deg + x, rad_lat, rad_eps);
            x - frac * (90.0 + ad)
        }
    };
    let (mut lo_x, mut hi_x) = (0.0, hi);
    if g(lo_x) > 0.0 {
        return None;
    }
    for _ in 0..PLACIDUS_BISECTION_ITERATIONS {
        let mid = 0.5 * (lo_x + hi_x);
        if g(mid) <= 0.0 {
            lo_x = mid;
        } else {
            hi_x = mid;
        }
    }
    let x = 0.5 * (lo_x + hi_x);
    let ra = if from_ic {
        wrap360(armc_deg + 180.0 - x)
    } else {
        wrap360(armc_deg + x)
    };
    let ra_rad = ra.to_radians();
    let (sin_ra, cos_ra) = sin_cos(ra_rad);
    let y = sin_ra;
    let xx = cos_ra * rad_eps.cos();
    Some(wrap360(y.atan2(xx).to_degrees()))
}

/// Placidus houses: the most
/// popular modern-Western system — time-based trisection of the
/// horizon→meridian travel (de Titis, *Primum Mobile*, 1657; Holden
/// 1977).
///
/// With `AD = arcsin(tan φ·tan δ)` and the semi-arc `SA = 90° ± AD`:
/// `H₁₁ = (90 + AD)/3`, `H₁₂ = 2(90 + AD)/3`,
/// `H₂ = 2(90 − AD)/3`, `H₃ = (90 − AD)/3`;
/// `RA₁₁/₁₂ = ARMC + H`, `RA₂/₃ = ARMC + 180° − H`. The iteration
/// recomputes `tan δ = sin RA·tan ε`, then `AD`, then the next `RA`
/// until `|ΔRA| < 1e-7°`, and each converged right ascension maps to
/// ecliptic longitude `λ = atan2(sin RA, cos RA·cos ε)`. Fast
/// fixed-point solves first, the bisection fallback where they diverge
/// within float noise of the polar circle, and the internal Porphyry
/// safety net when even the bisection bracket fails (unreachable with
/// the clamped AD).
///
/// All angles are degrees: `armc_deg` is the active RAMC (the entry
/// prologue owns the raising-selector polar-circle refusal and the
/// MC-flip; this kernel takes it as given), `asc_deg`/`mc_deg` seed
/// houses 1/10. Near the polar circle (threshold ≈ 90° − ε) the diurnal
/// arcs go circumpolar and the construction fails — the entry refuses
/// with `PolarCircle` and the fallback entry serves Porphyry. Returns
/// `[f64; 12]` houses 1–12.
///
/// Per-cusp `degnorm` shaping is the separate [`degnorm_deg`] step at
/// the public entry (see the module contract above).
pub fn houses_placidus(
    armc_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
    asc_deg: f64,
    mc_deg: f64,
) -> [f64; 12] {
    let rad_lat = lat_deg.to_radians();
    let rad_eps = eps_deg.to_radians();
    let c11 = iterate_placidus(armc_deg, rad_lat, rad_eps, 30.0)
        .or_else(|| bisect_placidus(armc_deg, rad_lat, rad_eps, 30.0));
    let c12 = iterate_placidus(armc_deg, rad_lat, rad_eps, 60.0)
        .or_else(|| bisect_placidus(armc_deg, rad_lat, rad_eps, 60.0));
    let c2 = iterate_placidus(armc_deg, rad_lat, rad_eps, 120.0)
        .or_else(|| bisect_placidus(armc_deg, rad_lat, rad_eps, 120.0));
    let c3 = iterate_placidus(armc_deg, rad_lat, rad_eps, 150.0)
        .or_else(|| bisect_placidus(armc_deg, rad_lat, rad_eps, 150.0));
    match (c11, c12, c2, c3) {
        (Some(c11), Some(c12), Some(c2), Some(c3)) => {
            let mut cusps = init_cardinal_cusps_13(asc_deg, mc_deg);
            cusps[11] = c11;
            cusps[12] = c12;
            cusps[2] = c2;
            cusps[3] = c3;
            cusps[5] = wrap360(c11 + 180.0);
            cusps[6] = wrap360(c12 + 180.0);
            cusps[8] = wrap360(c2 + 180.0);
            cusps[9] = wrap360(c3 + 180.0);
            strip_index0(cusps)
        }
        _ => {
            let porphyry = houses_porphyry(asc_deg, mc_deg);
            let mut full = [0.0_f64; 13];
            full[1..13].copy_from_slice(&porphyry);
            strip_index0(full)
        }
    }
}

/// Quadrant/time-based dispatch over the quadrant/time-based cusp-layer
/// systems with the release `houses`/`houses_armc` branch order and the
/// unknown-selector Placidus default (`else` arm, position
/// after the owned codes).
///
/// `hsys` is folded here ([`fold_hsys_case`]), like the release folds
/// before dispatching, so lowercase selectors (`'r'`, `'k'`, `'p'`, …)
/// resolve exactly as their uppercase forms. Koch (`K`, the served
/// Koch-kernel arm) is
/// served here with its kernel-level Porphyry fallback; every
/// other code outside the seven owned systems — the simple-ecliptic-layer
/// systems, the topocentric/specialized-layer systems (served by
/// [`houses_specialized`]), Gauquelin `G`/`g` (Gauquelin-sector layer; release
/// computes the 36 sectors for both, returning the first 12 for `g`,
/// so neither may take the Placidus default) — fails with
/// [`Error::HouseSystemUnserved`]
/// instead of silently returning Placidus cusps. Only truly unknown
/// selectors reach the default. The full `houses_armc` entry (ARMC-entry
/// layer) owns the complete dispatch. The polar-circle refusal mirrors the
/// entry prologue: raising selectors ([`polar_raise_applies`]) inside
/// the polar circle ([`is_polar_circle`]) fail with
/// [`Error::PolarCircle`] carrying
/// the `90 − eps` threshold, never with quiet fallback cusps.
pub fn houses_quadrant(
    armc_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
    asc_deg: f64,
    mc_deg: f64,
    hsys: char,
) -> Result<[f64; 12]> {
    let code = fold_hsys_case(hsys);
    if polar_raise_applies(code) && is_polar_circle(lat_deg, eps_deg) {
        return Err(Error::PolarCircle {
            latitude: lat_deg,
            threshold: polar_latitude_threshold_deg(eps_deg),
            obliquity: eps_deg,
            house_system: code,
        });
    }
    Ok(match code {
        'R' => houses_regiomontanus(armc_deg, lat_deg, eps_deg, asc_deg, mc_deg),
        'C' => houses_campanus(armc_deg, lat_deg, eps_deg, asc_deg, mc_deg),
        'B' => houses_alcabitius(armc_deg, lat_deg, eps_deg, asc_deg, mc_deg),
        'M' => houses_morinus(armc_deg, lat_deg, eps_deg, asc_deg, mc_deg),
        'X' => houses_meridian(armc_deg, lat_deg, eps_deg, asc_deg, mc_deg),
        'P' => houses_placidus(armc_deg, lat_deg, eps_deg, asc_deg, mc_deg),
        'K' => houses_koch(armc_deg, lat_deg, eps_deg, asc_deg, mc_deg),
        other if other == 'g' || KNOWN_HSYS_CODES.contains(other) => {
            return Err(Error::HouseSystemUnserved {
                house_system: other,
            });
        }
        _ => houses_placidus(armc_deg, lat_deg, eps_deg, asc_deg, mc_deg),
    })
}

/// Unit vector of a direction given by right ascension/declination
/// (`_unit_radec` ). Both same-argument pairs go through
/// [`sin_cos`] (separate-evaluation guard: the
/// declination pair fused into a 1-ulp-off platform `sincos` in release
/// builds, demonstrated by two APC sweep cases; the right-ascension
/// pair was already guarded).
pub(super) fn unit_radec(ra_rad: f64, dec_rad: f64) -> EqVec {
    let (sin_dec, cos_dec) = sin_cos(dec_rad);
    let (sin_ra, cos_ra) = sin_cos(ra_rad);
    [cos_dec * cos_ra, cos_dec * sin_ra, sin_dec]
}

/// Right-handed cross product (`_cross` ).
pub(super) fn vec_cross(u: EqVec, v: EqVec) -> EqVec {
    [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ]
}

/// Euclidean dot product (`_dot` ).
pub(super) fn vec_dot(u: EqVec, v: EqVec) -> f64 {
    u[0] * v[0] + u[1] * v[1] + u[2] * v[2]
}

/// Unit-length rescaling (`_normalized`; callers exclude
/// degenerate input, as in the release).
pub(super) fn vec_normalized(v: EqVec) -> EqVec {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    [v[0] / n, v[1] / n, v[2] / n]
}

/// Ecliptic longitude of an equatorial direction, degrees in `[0, 360)`
/// (`_ecliptic_lon_deg` ).
pub(super) fn ecliptic_lon_deg(v: EqVec, eps_rad: f64) -> f64 {
    let (sin_eps, cos_eps) = sin_cos(eps_rad);
    let lon = (v[1] * cos_eps + v[2] * sin_eps).atan2(v[0]);
    wrap360(lon.to_degrees())
}

/// Longitude where the hour circle through `v` meets the ecliptic
/// (`_hour_circle_lon_deg` ): the hour circle runs through both
/// celestial poles and `v`, crossing the ecliptic at `v`'s right
/// ascension.
pub(super) fn hour_circle_lon_deg(v: EqVec, eps_rad: f64) -> f64 {
    let ra = v[1].atan2(v[0]);
    let (sin_ra, cos_ra) = sin_cos(ra);
    let lon = sin_ra.atan2(cos_ra * eps_rad.cos());
    wrap360(lon.to_degrees())
}

/// Zenith, east point and north point of the horizon as unit vectors
/// (`_horizon_frame` ).
pub(super) fn horizon_frame(armc_rad: f64, lat_rad: f64) -> (EqVec, EqVec, EqVec) {
    let zenith = unit_radec(armc_rad, lat_rad);
    let east = unit_radec(armc_rad + std::f64::consts::FRAC_PI_2, 0.0);
    let north = vec_cross(zenith, east);
    (zenith, east, north)
}

/// Ecliptic longitude cut by the great circle through `axis` and
/// `anchor`.
///
/// `axis` is the unit vector of the horizon's north point (the circle
/// also contains the south point `-axis`). Of the two ecliptic
/// intersections `same_side` selects the one on the anchor's half of
/// the circle; otherwise the antipodal half. Two degeneracy fallbacks:
/// an axis-parallel circle carries the anchor's hour circle, and an
/// in-ecliptic circle reports the anchor's own longitude (each with the
/// `+180°` antipode when `same_side` is false). Branch selection
/// compares the `dot(point, half)` signs — a geometric half-circle
/// invariant, stable across the 0/360 seam, both hemispheres and polar
/// frames — never longitude quadrants.
///
/// `anchor`/`axis` are equatorial unit vectors, `eps_rad` radians;
/// returns degrees in `[0, 360)`; scalar on both sides.
pub(super) fn house_circle_lon_deg(
    anchor: EqVec,
    axis: EqVec,
    eps_rad: f64,
    same_side: bool,
) -> f64 {
    let normal = vec_cross(axis, anchor);
    if vec_dot(normal, normal) < HOUSE_CIRCLE_TINY {
        let longitude = hour_circle_lon_deg(anchor, eps_rad);
        return if same_side {
            longitude
        } else {
            wrap360(longitude + 180.0)
        };
    }
    let ecl_pole: EqVec = {
        let (sin_eps, cos_eps) = sin_cos(eps_rad);
        [0.0, -sin_eps, cos_eps]
    };
    let mut meet = vec_cross(normal, ecl_pole);
    if vec_dot(meet, meet) < HOUSE_CIRCLE_TINY {
        let longitude = ecliptic_lon_deg(anchor, eps_rad);
        return if same_side {
            longitude
        } else {
            wrap360(longitude + 180.0)
        };
    }
    let half = vec_cross(normal, axis);
    let side_product = vec_dot(meet, half) * vec_dot(anchor, half);
    if (side_product < 0.0) == same_side {
        meet = [-meet[0], -meet[1], -meet[2]];
    }
    ecliptic_lon_deg(meet, eps_rad)
}

/// Hour angles of the eight intermediate division points of a declination
/// parallel (`_semiarc_cusp_hour_angles` ).
///
/// `day_semiarc_rad` is half the diurnal arc. Cusps 12/11 and 9/8
/// trisect the diurnal semi-arcs east and west of the meridian; cusps
/// 2/3 and 5/6 trisect the nocturnal semi-arcs on either side of the
/// lower meridian. Hour angle is west-positive.
pub(super) fn semiarc_cusp_hour_angles(day_semiarc_rad: f64) -> [(usize, f64); 8] {
    let day = day_semiarc_rad;
    let night = std::f64::consts::PI - day;
    [
        (11, -day / 3.0),
        (12, -2.0 * day / 3.0),
        (2, -day - night / 3.0),
        (3, -day - 2.0 * night / 3.0),
        (9, day / 3.0),
        (8, 2.0 * day / 3.0),
        (6, day + night / 3.0),
        (5, day + 2.0 * night / 3.0),
    ]
}

/// Signed angular distance `p1 − p2` normalized to `[-180, 180]`
/// e.g.
/// `(10, 20) → −10`, `(350, 10) → −20`, `(10, 350) → 20`,
/// `(180, 0) → −180`, via the `% 360` fold with the `>= 180 → −360`
/// rule.
///
/// Note the `>= 180` fold: an exact half-turn reports `-180.0`, unlike
/// [`angular_diff`] (strict `> 180`, reporting `+180.0`). Inventoried with
/// the system-inventory layer, implemented here where the APC kernel needs
/// it. Scalar on both sides.
pub fn difdeg2n_deg(p1: f64, p2: f64) -> f64 {
    crate::utils::difdeg2n(p1, p2)
}

/// Ecliptic longitude from ARMC and pole height (Smart/Meeus rising-point construction): computes
/// the obliquity sine/cosine internally and delegates to
/// [`ra_to_ecliptic_longitude`].
///
/// `armc_deg` and `pole_height_deg` are degrees in, ecliptic longitude
/// in `[0, 360)` out; `eps_deg` is degrees; `lat_deg` is unused, kept
/// for call-site compatibility like the release signature. Shared
/// ARMC-entry ASCMC primitive (also the Koch/Sunshine rising-point
/// construction). Scalar on both sides.
pub fn calc_ascendant(armc_deg: f64, eps_deg: f64, _lat_deg: f64, pole_height_deg: f64) -> f64 {
    let eps_rad = eps_deg.to_radians();
    let (sin_eps, cos_eps) = sin_cos(eps_rad);
    ra_to_ecliptic_longitude(armc_deg, pole_height_deg, sin_eps, cos_eps)
}

/// Polich/Page topocentric houses: the topocentric system with the observer-surface perspective
/// (Polich & Page, *Spica*, 1964) — Regiomontanus-like modified poles
/// with `tan(P) = tan(φ)·factor` (`1/3` → houses 11/3, `2/3` → houses
/// 12/2), `R = ARMC + H − 90°` and
/// `λ = atan2(cos R, −(sin R·cos ε + tan P·sin ε))` for
/// `H = 30°/60°/120°/150°` (houses 11/12/2/3); 5/6/8/9 are the 180°
/// opposites.
///
/// Arguments are degrees (active ARMC, latitude, true obliquity,
/// precomputed Ascendant/MC seeding houses 1/10); takes the passed
/// angles as-is (entry owns the MC-flip and ASCMC assembly). Returns
/// `[f64; 12]` houses 1–12.
pub fn houses_polich_page(
    armc_deg: f64,
    lat_deg: f64,
    eps_deg: f64,
    asc_deg: f64,
    mc_deg: f64,
) -> [f64; 12] {
    fn calc_cusp(armc_deg: f64, rad_lat: f64, rad_eps: f64, offset_deg: f64, factor: f64) -> f64 {
        let tan_pole = rad_lat.tan() * factor;
        let r_deg = wrap360(armc_deg + offset_deg - 90.0);
        let r_rad = r_deg.to_radians();
        let (sin_r, cos_r) = sin_cos(r_rad);
        let (sin_eps, cos_eps) = sin_cos(rad_eps);
        let num = cos_r;
        let den = -(sin_r * cos_eps + tan_pole * sin_eps);
        wrap360(num.atan2(den).to_degrees())
    }

    let mut cusps = init_cardinal_cusps_13(asc_deg, mc_deg);
    let rad_lat = lat_deg.to_radians();
    let rad_eps = eps_deg.to_radians();
    cusps[11] = calc_cusp(armc_deg, rad_lat, rad_eps, 30.0, 1.0 / 3.0);
    cusps[12] = calc_cusp(armc_deg, rad_lat, rad_eps, 60.0, 2.0 / 3.0);
    cusps[2] = calc_cusp(armc_deg, rad_lat, rad_eps, 120.0, 2.0 / 3.0);
    cusps[3] = calc_cusp(armc_deg, rad_lat, rad_eps, 150.0, 1.0 / 3.0);
    set_opposite_cusps_13(&mut cusps);
    strip_index0(cusps)
}

/// Carter poli-equatorial houses:
/// equal 30° equatorial divisions from the Ascendant's right ascension
/// with hour-circle projection (Carter, *Essays on the Foundations of
/// Astrology*): `RA_Asc = atan2(sin λ_Asc·cos ε, cos λ_Asc)`,
/// `RA = RA_Asc + (i−1)·30°`,
/// `λ = atan2(sin RA, cos RA·cos ε)`.
///
/// Arguments are degrees (true obliquity, Ascendant; ARMC/latitude/MC
/// unused — the call shape keeps them); takes the passed angles as-is.
/// Returns `[f64; 12]` houses 1–12.
pub fn houses_carter(
    _armc_deg: f64,
    _lat_deg: f64,
    eps_deg: f64,
    asc_deg: f64,
    _mc_deg: f64,
) -> [f64; 12] {
    fn lon_from_ra(rad_eps: f64, ra_deg: f64) -> f64 {
        let ra_rad = ra_deg.to_radians();
        let (sin_ra, cos_ra) = sin_cos(ra_rad);
        let y = sin_ra;
        let x = cos_ra * rad_eps.cos();
        wrap360(y.atan2(x).to_degrees())
    }

    let rad_eps = eps_deg.to_radians();
    let asc_rad = asc_deg.to_radians();
    let (sin_asc, cos_asc) = sin_cos(asc_rad);
    let y = sin_asc * rad_eps.cos();
    let x = cos_asc;
    let ra_asc = wrap360(y.atan2(x).to_degrees());
    let mut cusps = [0.0_f64; 13];
    for (n, cusp) in cusps.iter_mut().enumerate().skip(1) {
        let ra = wrap360(ra_asc + (n as f64 - 1.0) * 30.0);
        *cusp = lon_from_ra(rad_eps, ra);
    }
    strip_index0(cusps)
}
