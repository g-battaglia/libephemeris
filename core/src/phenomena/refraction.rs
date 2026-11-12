//! Atmospheric refraction and horizon dip.
//!
//! The ICAO atmosphere ray tracer evaluates the Green (1985, eq. 4.25)
//! refraction integral with 120-point Gauss–Legendre quadrature. A compact
//! Sæmundsson/Meeus model supplies the public altitude correction, and the
//! Thom model supplies the refraction-corrected horizon dip.
//!
//! Altitudes and corrections are in degrees; pressure is in mbar, heights
//! and radii are in metres, and each entry identifies Celsius or Kelvin.
//! The functions require no ephemeris or clock. Pressure estimation from
//! observer altitude is an explicit operation, not an implicit kernel step.
//!
//! [`AtmosphereProfile`] owns the sampled atmosphere and can be reused by
//! its caller. Pure scalar kernels return IEEE binary64 values; their
//! non-physical-input branches and iteration limits are documented per entry.

use crate::refraction_tables::{QUAD_N, QUAD_NODE_BITS, QUAD_WEIGHT_BITS};

/// Direction flag: convert true altitude to apparent altitude
pub const TRUE_TO_APP: i32 = 0;
/// Direction flag: convert apparent altitude to true altitude
pub const APP_TO_TRUE: i32 = 1;

/// Mean Earth radius in metres.
pub const EARTH_RADIUS_M: f64 = 6_371_000.0;
/// Standard gravity in m/s².
pub const STANDARD_GRAVITY_M_S2: f64 = 9.80665;
/// Molar mass of dry air in kg/mol.
pub const MOLAR_MASS_AIR_KG_MOL: f64 = 0.0289644;
/// Universal gas constant in J/(mol·K).
pub const GAS_CONSTANT_J_MOL_K: f64 = 8.314462618;
/// Hydrostatic exponent `g·M/R` in K/m.
pub const HYDROSTATIC_EXPONENT_K_M: f64 =
    STANDARD_GRAVITY_M_S2 * MOLAR_MASS_AIR_KG_MOL / GAS_CONSTANT_J_MOL_K;
/// Refractive-index coefficient.
pub const REFRACTIVE_INDEX_COEFF: f64 = 8.060e-5;
/// Top of the modeled atmosphere in metres.
pub const ATMOSPHERE_TOP_M: f64 = 84_852.0;
/// Standard sea-level pressure in mbar.
pub const STANDARD_PRESSURE_MBAR: f64 = 1013.25;
/// Standard sea-level temperature in K.
pub const STANDARD_TEMPERATURE_K: f64 = 288.15;

/// ICAO atmosphere layers as `(base_altitude_m, lapse_rate_K_per_m,
/// T_base_K)`.
pub const ICAO_LAYERS: [(f64, f64, f64); 7] = [
    (0.0, -0.0065, 288.15),
    (11_000.0, 0.0, 216.65),
    (20_000.0, 0.001, 216.65),
    (32_000.0, 0.0028, 228.65),
    (47_000.0, 0.0, 270.65),
    (51_000.0, -0.0028, 270.65),
    (71_000.0, -0.002, 214.65),
];

/// Conventional compat pressure reference in mbar
pub const COMPAT_PRESSURE_REF_MBAR: f64 = 1010.0;
/// Conventional compat temperature-scale numerator in K
pub const COMPAT_TEMPERATURE_NUMERATOR_K: f64 = 283.0;
/// Conventional compat temperature-scale offset in K
pub const COMPAT_TEMPERATURE_OFFSET_K: f64 = 273.0;
/// Compat forward low-altitude cutoff in degrees
pub const COMPAT_TRUE_LOW_LIMIT_DEG: f64 = -5.0;
/// Compat forward branch altitude in degrees
pub const COMPAT_TRUE_BRANCH_DEG: f64 = 15.0;
/// Compat reverse pole offset in apparent-altitude degrees: the Bennett
/// closed form divides by `apparent_alt + 4.4`, and an exactly-zero
/// divisor returns the input unchanged.
pub const COMPAT_REVERSE_POLE_OFFSET_DEG: f64 = 4.4;
/// Compat reverse pole numerator in degrees.
pub const COMPAT_REVERSE_POLE_NUMERATOR_DEG: f64 = 7.31;
/// Compat reverse zenith guard in degrees.
pub const COMPAT_REVERSE_ZENITH_LIMIT_DEG: f64 = 90.0 - 1e-10;
/// Newton iteration budget for the ray-traced apparent-to-true inversion
pub const INVERSE_NEWTON_ITERATIONS: usize = 8;
/// Newton residual stop band in degrees.
pub const INVERSE_NEWTON_TOLERANCE_DEG: f64 = 1e-12;
/// Newton forward-difference step in degrees.
pub const INVERSE_NEWTON_STEP_DEG: f64 = 0.001;
/// Newton derivative floor in apparent-degrees-per-true-degree: iteration
/// stops instead of dividing by a smaller slope.
pub const INVERSE_NEWTON_DERIVATIVE_FLOOR: f64 = 1e-15;

/// Equatorial Earth radius in metres used by the geometric dip
pub const DIP_EARTH_RADIUS_M: f64 = 6_378_136.6;

/// Below-horizon linear-extrapolation pivot in true-zenith degrees
pub const TRACE_EXTRAPOLATION_PIVOT_DEG: f64 = 91.0;

/// Fixed-point iteration budget for the `z_obs = z_true - R(z_obs)` loop
pub const TRACE_ITERATIONS: usize = 8;
/// Fixed-point stop band in degrees.
pub const TRACE_TOLERANCE_DEG: f64 = 1e-12;
/// Central-difference half-step in metres for `dn/dr`
pub const DN_DR_STEP_M: f64 = 1.0;

/// Temperature at altitude `h_m` in metres per the ICAO profile, in Kelvin
pub fn icao_temperature_k(h_m: f64) -> f64 {
    let mut i = ICAO_LAYERS.len();
    while i > 0 {
        i -= 1;
        let (h_base, lapse, t_base) = ICAO_LAYERS[i];
        if h_m >= h_base {
            return t_base + lapse * (h_m - h_base);
        }
    }
    ICAO_LAYERS[0].2
}

/// Pressure at altitude `h_m` in metres per the ICAO profile, in mbar
pub fn icao_pressure_mbar(h_m: f64) -> f64 {
    let mut p = STANDARD_PRESSURE_MBAR;
    let mut i = 0;
    while i < ICAO_LAYERS.len() {
        let (h_base, lapse, t_base) = ICAO_LAYERS[i];
        let h_top = if i + 1 < ICAO_LAYERS.len() {
            ICAO_LAYERS[i + 1].0
        } else {
            ATMOSPHERE_TOP_M
        };
        if h_m <= h_top {
            let dh = h_m - h_base;
            if lapse.abs() < 1e-10 {
                p *= (-HYDROSTATIC_EXPONENT_K_M * dh / t_base).exp();
            } else {
                p *= ((t_base + lapse * dh) / t_base).powf(-HYDROSTATIC_EXPONENT_K_M / lapse);
            }
            return p;
        }
        let dh = h_top - h_base;
        if lapse.abs() < 1e-10 {
            p *= (-HYDROSTATIC_EXPONENT_K_M * dh / t_base).exp();
        } else {
            p *= ((t_base + lapse * dh) / t_base).powf(-HYDROSTATIC_EXPONENT_K_M / lapse);
        }
        i += 1;
    }
    p
}

/// Refractive index of air at `p_mbar` and `t_kelvin`
/// Non-physical inputs (non-positive pressure or temperature) report no refractivity, i.e. exactly 1.0.
pub fn refractive_index(p_mbar: f64, t_kelvin: f64) -> f64 {
    if t_kelvin <= 0.0 || p_mbar <= 0.0 {
        return 1.0;
    }
    1.0 + REFRACTIVE_INDEX_COEFF * p_mbar / t_kelvin
}

/// Temperature at altitude `h_m` under the observer's surface conditions
/// in the troposphere and the shifted ICAO profile above 11 km, in Kelvin
pub fn modified_temperature_k(
    h_m: f64,
    obs_alt_m: f64,
    surface_t_k: f64,
    lapse_rate_k_m: f64,
) -> f64 {
    let tropo_top = 11_000.0;
    if h_m < tropo_top {
        return surface_t_k - lapse_rate_k_m * (h_m - obs_alt_m);
    }
    let t_at_tropo = surface_t_k - lapse_rate_k_m * (tropo_top - obs_alt_m);
    let mut i = 1;
    while i < ICAO_LAYERS.len() {
        let h_base = ICAO_LAYERS[i].0;
        let lapse_i = ICAO_LAYERS[i].1;
        let h_top = if i + 1 < ICAO_LAYERS.len() {
            ICAO_LAYERS[i + 1].0
        } else {
            ATMOSPHERE_TOP_M
        };
        if h_m <= h_top {
            if i == 1 {
                return t_at_tropo + lapse_i * (h_m - h_base);
            }
            let offset = t_at_tropo - ICAO_LAYERS[1].2;
            return ICAO_LAYERS[i].2 + offset + lapse_i * (h_m - h_base);
        }
        i += 1;
    }
    200.0
}

/// Pressure at altitude `h_m` under the modified atmosphere, in mbar
pub fn modified_pressure_mbar(
    h_m: f64,
    obs_alt_m: f64,
    obs_p_mbar: f64,
    surface_t_k: f64,
    lapse_rate_k_m: f64,
) -> f64 {
    let tropo_top = 11_000.0;
    if h_m <= tropo_top {
        let t_obs = surface_t_k;
        let t_h = t_obs - lapse_rate_k_m * (h_m - obs_alt_m);
        if lapse_rate_k_m.abs() < 1e-10 {
            return obs_p_mbar * (-HYDROSTATIC_EXPONENT_K_M * (h_m - obs_alt_m) / t_obs).exp();
        }
        return obs_p_mbar * (t_h / t_obs).powf(HYDROSTATIC_EXPONENT_K_M / lapse_rate_k_m);
    }

    let t_obs = surface_t_k;
    let t_tropo = t_obs - lapse_rate_k_m * (tropo_top - obs_alt_m);
    let p_tropo = if lapse_rate_k_m.abs() < 1e-10 {
        obs_p_mbar * (-HYDROSTATIC_EXPONENT_K_M * (tropo_top - obs_alt_m) / t_obs).exp()
    } else {
        obs_p_mbar * (t_tropo / t_obs).powf(HYDROSTATIC_EXPONENT_K_M / lapse_rate_k_m)
    };

    let mut p = p_tropo;
    let t_at_tropo = t_tropo;
    let mut i = 1;
    while i < ICAO_LAYERS.len() {
        let h_base = ICAO_LAYERS[i].0;
        let lapse_i = ICAO_LAYERS[i].1;
        let h_top_i = if i + 1 < ICAO_LAYERS.len() {
            ICAO_LAYERS[i + 1].0
        } else {
            ATMOSPHERE_TOP_M
        };
        let t_base_i = if i == 1 {
            t_at_tropo
        } else {
            ICAO_LAYERS[i].2 + (t_at_tropo - ICAO_LAYERS[1].2)
        };
        if h_m <= h_top_i {
            let dh = h_m - h_base;
            if lapse_i.abs() < 1e-10 {
                p *= (-HYDROSTATIC_EXPONENT_K_M * dh / t_base_i).exp();
            } else {
                p *= ((t_base_i + lapse_i * dh) / t_base_i)
                    .powf(-HYDROSTATIC_EXPONENT_K_M / lapse_i);
            }
            return p;
        }
        let dh = h_top_i - h_base;
        if lapse_i.abs() < 1e-10 {
            p *= (-HYDROSTATIC_EXPONENT_K_M * dh / t_base_i).exp();
        } else {
            p *= ((t_base_i + lapse_i * dh) / t_base_i).powf(-HYDROSTATIC_EXPONENT_K_M / lapse_i);
        }
        i += 1;
    }
    p
}

/// Refractive index at altitude `h_m` in the modified atmosphere
pub fn n_at_height(
    h_m: f64,
    obs_alt_m: f64,
    obs_p_mbar: f64,
    obs_t_k: f64,
    lapse_rate_k_m: f64,
) -> f64 {
    let t = modified_temperature_k(h_m, obs_alt_m, obs_t_k, lapse_rate_k_m);
    let p = modified_pressure_mbar(h_m, obs_alt_m, obs_p_mbar, obs_t_k, lapse_rate_k_m);
    refractive_index(p, t)
}

/// Derivative `dn/dr` at altitude `h_m` via the ±1 m central difference
pub fn dn_dr_at_height(
    h_m: f64,
    obs_alt_m: f64,
    obs_p_mbar: f64,
    obs_t_k: f64,
    lapse_rate_k_m: f64,
) -> f64 {
    let n_plus = n_at_height(
        h_m + DN_DR_STEP_M,
        obs_alt_m,
        obs_p_mbar,
        obs_t_k,
        lapse_rate_k_m,
    );
    let n_minus = n_at_height(
        h_m - DN_DR_STEP_M,
        obs_alt_m,
        obs_p_mbar,
        obs_t_k,
        lapse_rate_k_m,
    );
    (n_plus - n_minus) / (2.0 * DN_DR_STEP_M)
}

/// Precomputed atmosphere profile at the Gauss-Legendre quadrature points
///
/// Fields are private: the profile is built by [`AtmosphereProfile::new`]
/// and consumed by [`refraction_integral`]. Unlike the release
/// single-slot global cache, profiles are explicit values with no shared
/// mutable state.
/// Radius maps as r = r_obs + L*u^2 with L = r_top - r_obs, so dr = 2*L*u*du desingularizes the horizon integrand; every quadrature point stores its radius, index, and gradient up front.
pub struct AtmosphereProfile {
    n_obs: f64,
    r_obs: f64,
    l_m: f64,
    q_r: [f64; QUAD_N],
    q_n: [f64; QUAD_N],
    q_dn_dr: [f64; QUAD_N],
}

impl AtmosphereProfile {
    /// Build the profile for the observer conditions.
    pub fn new(obs_alt_m: f64, obs_p_mbar: f64, obs_t_k: f64, lapse_rate_k_m: f64) -> Self {
        let r_obs = EARTH_RADIUS_M + obs_alt_m;
        let r_top = EARTH_RADIUS_M + ATMOSPHERE_TOP_M;
        let l_m = r_top - r_obs;
        let mut q_r = [0.0f64; QUAD_N];
        let mut q_n = [0.0f64; QUAD_N];
        let mut q_dn_dr = [0.0f64; QUAD_N];
        let mut k = 0;
        while k < QUAD_N {
            let node = f64::from_bits(QUAD_NODE_BITS[k]);
            let u = 0.5 * (node + 1.0);
            let r = r_obs + l_m * u * u;
            let h = r - EARTH_RADIUS_M;
            q_r[k] = r;
            q_n[k] = n_at_height(h, obs_alt_m, obs_p_mbar, obs_t_k, lapse_rate_k_m);
            q_dn_dr[k] = dn_dr_at_height(h, obs_alt_m, obs_p_mbar, obs_t_k, lapse_rate_k_m);
            k += 1;
        }
        let n_obs = n_at_height(obs_alt_m, obs_alt_m, obs_p_mbar, obs_t_k, lapse_rate_k_m);
        Self {
            n_obs,
            r_obs,
            l_m,
            q_r,
            q_n,
            q_dn_dr,
        }
    }

    /// Refractive index at the observer.
    pub fn n_obs(&self) -> f64 {
        self.n_obs
    }

    /// Geocentric observer radius in metres.
    pub fn r_obs(&self) -> f64 {
        self.r_obs
    }
}

/// Evaluate the refraction integral for a ray arriving at the observer at
/// apparent zenith angle `observed_zenith_deg`; returns refraction in
/// degrees `>= 0`.
/// The Bouguer invariant C = n_obs*r_obs*sin(z_obs) sets each term; quadrature points the ray never reaches (non-positive argument) are skipped, and the `[-1,1]`-to-`[0,1]` mapping contributes the one-half weight.
pub fn refraction_integral(observed_zenith_deg: f64, profile: &AtmosphereProfile) -> f64 {
    let z_rad = observed_zenith_deg.to_radians();
    let sin_z = z_rad.sin();
    let c = profile.n_obs * profile.r_obs * sin_z;
    let c2 = c * c;
    let l_m = profile.l_m;

    let mut total = 0.0;
    let mut k = 0;
    while k < QUAD_N {
        let node = f64::from_bits(QUAD_NODE_BITS[k]);
        let weight = f64::from_bits(QUAD_WEIGHT_BITS[k]);
        let u = 0.5 * (node + 1.0);
        let r = profile.q_r[k];
        let n = profile.q_n[k];
        let dn_dr = profile.q_dn_dr[k];

        let nr = n * r;
        let arg = nr * nr / c2 - 1.0;
        if arg <= 0.0 {
            k += 1;
            continue;
        }
        let f = -dn_dr / (n * arg.sqrt()) * 2.0 * l_m * u;
        total += 0.5 * weight * f;
        k += 1;
    }

    (0.0f64).max(total.to_degrees())
}

/// Refraction for a source at true zenith angle `true_zenith_deg`, in
/// degrees `>= 0`.
/// Non-positive zenith, non-positive pressure, and non-physical observer temperatures report zero; past 91 deg true zenith the 90/91 deg pair sets a linear extension, and the observed-zenith fixed point runs at most 8 steps to a 1e-12 deg stop, clamped at zero.
#[allow(clippy::too_many_arguments)]
pub fn trace_ray(
    true_zenith_deg: f64,
    obs_alt_m: f64,
    obs_pressure_mbar: f64,
    obs_temperature_c: f64,
    lapse_rate_k_m: f64,
) -> f64 {
    if true_zenith_deg <= 0.0 {
        return 0.0;
    }
    if obs_pressure_mbar <= 0.0 {
        return 0.0;
    }
    let obs_t_k = obs_temperature_c + 273.15;
    if obs_t_k <= 0.0 {
        return 0.0;
    }
    if true_zenith_deg > TRACE_EXTRAPOLATION_PIVOT_DEG {
        let r_91 = trace_ray(
            TRACE_EXTRAPOLATION_PIVOT_DEG,
            obs_alt_m,
            obs_pressure_mbar,
            obs_temperature_c,
            lapse_rate_k_m,
        );
        let r_90 = trace_ray(
            90.0,
            obs_alt_m,
            obs_pressure_mbar,
            obs_temperature_c,
            lapse_rate_k_m,
        );
        let slope = r_91 - r_90;
        return (0.0f64).max(r_91 + slope * (true_zenith_deg - TRACE_EXTRAPOLATION_PIVOT_DEG));
    }

    let profile = AtmosphereProfile::new(obs_alt_m, obs_pressure_mbar, obs_t_k, lapse_rate_k_m);
    let mut z_obs = true_zenith_deg;
    let mut r = 0.0;
    let mut step = 0;
    while step < TRACE_ITERATIONS {
        r = refraction_integral(z_obs, &profile);
        let z_obs_new = true_zenith_deg - r;
        if (z_obs_new - z_obs).abs() < TRACE_TOLERANCE_DEG {
            break;
        }
        z_obs = z_obs_new;
        step += 1;
    }
    (0.0f64).max(r)
}

/// ICAO ray-traced refraction for a true (geometric) altitude in degrees;
/// returns refraction in degrees `>= 0`
pub fn calc_refraction_true_to_app(
    true_alt_deg: f64,
    pressure_mbar: f64,
    temperature_c: f64,
    obs_alt_m: f64,
    lapse_rate_k_m: f64,
) -> f64 {
    if pressure_mbar <= 0.0 {
        return 0.0;
    }
    let zenith = 90.0 - true_alt_deg;
    trace_ray(
        zenith,
        obs_alt_m,
        pressure_mbar,
        temperature_c,
        lapse_rate_k_m,
    )
}

/// ICAO ray-traced refraction for a true altitude; legacy internal alias
pub fn calc_refraction_ref_true_to_app(
    true_alt_deg: f64,
    pressure_mbar: f64,
    temperature_c: f64,
    obs_alt_m: f64,
    lapse_rate_k_m: f64,
) -> f64 {
    calc_refraction_true_to_app(
        true_alt_deg,
        pressure_mbar,
        temperature_c,
        obs_alt_m,
        lapse_rate_k_m,
    )
}

/// ICAO ray-traced refraction for an apparent (observed) altitude in
/// degrees; returns refraction in degrees `>= 0`.
///
/// Newton-Raphson inversion of the true-to-apparent mapping with the
/// same 8-step budget, the `0.001` forward-difference derivative,
/// the `1e-12` residual stop and the `1e-15` derivative floor. Like the
/// forward kernel, `pressure <= 0` reports no refraction and the result
/// is clamped at zero. No trigonometry here beyond what the forward
/// kernel evaluates (single sin/tan per site, never same-argument
/// pairs, so no platform sincos merging concern).
#[allow(clippy::too_many_arguments)]
pub fn calc_refraction_app_to_true(
    apparent_alt_deg: f64,
    pressure_mbar: f64,
    temperature_c: f64,
    obs_alt_m: f64,
    lapse_rate_k_m: f64,
) -> f64 {
    if pressure_mbar <= 0.0 {
        return 0.0;
    }
    let mut true_est = apparent_alt_deg;
    for _ in 0..INVERSE_NEWTON_ITERATIONS {
        let r = calc_refraction_true_to_app(
            true_est,
            pressure_mbar,
            temperature_c,
            obs_alt_m,
            lapse_rate_k_m,
        );
        let residual = true_est + r - apparent_alt_deg;
        if residual.abs() < INVERSE_NEWTON_TOLERANCE_DEG {
            break;
        }
        let r_plus = calc_refraction_true_to_app(
            true_est + INVERSE_NEWTON_STEP_DEG,
            pressure_mbar,
            temperature_c,
            obs_alt_m,
            lapse_rate_k_m,
        );
        let dapp_dtrue = 1.0 + (r_plus - r) / INVERSE_NEWTON_STEP_DEG;
        if dapp_dtrue.abs() < INVERSE_NEWTON_DERIVATIVE_FLOOR {
            break;
        }
        true_est -= residual / dapp_dtrue;
    }
    (0.0f64).max(apparent_alt_deg - true_est)
}

/// ICAO ray-traced refraction for an apparent altitude; legacy internal
/// alias.
#[allow(clippy::too_many_arguments)]
pub fn calc_refraction_ref_app_to_true(
    apparent_alt_deg: f64,
    pressure_mbar: f64,
    temperature_c: f64,
    obs_alt_m: f64,
    lapse_rate_k_m: f64,
) -> f64 {
    calc_refraction_app_to_true(
        apparent_alt_deg,
        pressure_mbar,
        temperature_c,
        obs_alt_m,
        lapse_rate_k_m,
    )
}

/// Conventional pressure/temperature refraction scale for the compat
/// closed forms.
/// A zero denominator never raises: a zero-or-NaN numerator yields NaN, otherwise signed infinity following the numerator.
pub fn compat_atmosphere_scale(pressure_mbar: f64, temperature_c: f64) -> f64 {
    let numerator = pressure_mbar / COMPAT_PRESSURE_REF_MBAR * COMPAT_TEMPERATURE_NUMERATOR_K;
    let denominator = COMPAT_TEMPERATURE_OFFSET_K + temperature_c;
    if denominator != 0.0 {
        return numerator / denominator;
    }
    if numerator == 0.0 || numerator.is_nan() {
        return f64::NAN;
    }
    f64::copysign(f64::INFINITY, numerator)
}

/// Plain-API apparent altitude for a true altitude under the Saemundsson
/// / high-altitude tangent-series closed forms. Inputs at or below −5 deg, and
/// candidates that do not rise strictly above zero, are returned
/// unchanged.
pub fn calc_refrac_compat_true_to_app(
    true_alt_deg: f64,
    pressure_mbar: f64,
    temperature_c: f64,
) -> f64 {
    if !true_alt_deg.is_finite() || true_alt_deg <= COMPAT_TRUE_LOW_LIMIT_DEG {
        return true_alt_deg;
    }
    let scale = compat_atmosphere_scale(pressure_mbar, temperature_c);
    let correction = if true_alt_deg <= COMPAT_TRUE_BRANCH_DEG {
        let angle = true_alt_deg + 10.3 / (true_alt_deg + 5.11);
        1.02 / angle.to_radians().tan() * scale / 60.0
    } else {
        let tangent = (90.0 - true_alt_deg).to_radians().tan();
        (58.276 * tangent - 0.0824 * tangent * tangent * tangent) * scale / 3600.0
    };
    let apparent = true_alt_deg + correction;
    if apparent > 0.0 {
        apparent
    } else {
        true_alt_deg
    }
}

/// Plain-API true altitude for an apparent altitude under Bennett's
/// closed form.
///
/// This is the independently published Bennett approximation (Meeus
/// 1998, ch. 16, eq. 16.3), not a numerical inversion of
/// [`calc_refrac_compat_true_to_app`]: the two directions are
/// independent fits and do not round-trip exactly. Non-finite inputs
/// pass through unchanged; the exact `-4.4` pole (a zero divisor that
/// the release reports via `ZeroDivisionError`), angles at or above the
/// reverse zenith limit, an exactly-zero tangent and candidates not
/// strictly above zero all return the input unchanged,.
pub fn calc_refrac_compat_app_to_true(
    apparent_alt_deg: f64,
    pressure_mbar: f64,
    temperature_c: f64,
) -> f64 {
    if !apparent_alt_deg.is_finite() {
        return apparent_alt_deg;
    }
    let denom = apparent_alt_deg + COMPAT_REVERSE_POLE_OFFSET_DEG;
    if denom == 0.0 {
        return apparent_alt_deg;
    }
    let angle = apparent_alt_deg + COMPAT_REVERSE_POLE_NUMERATOR_DEG / denom;
    if !angle.is_finite() || angle >= COMPAT_REVERSE_ZENITH_LIMIT_DEG {
        return apparent_alt_deg;
    }
    let tangent = angle.to_radians().tan();
    if tangent == 0.0 {
        return apparent_alt_deg;
    }
    let scale = compat_atmosphere_scale(pressure_mbar, temperature_c);
    let correction = 1.0 / tangent * scale / 60.0;
    let true_alt = apparent_alt_deg - correction;
    if true_alt > 0.0 {
        true_alt
    } else {
        apparent_alt_deg
    }
}

/// Thom refraction-corrected dip of the horizon for an observer at
/// `obs_alt_m` metres, in degrees (negative; release `calc_dip`).
/// Form: geometric acos(R/(R+h)) scaled by sqrt(max(0, 1 - 1.848*k*P/T^2)) with k = (0.0342 + lapse)/(0.154*0.0238) and T = 273.16 + attemp_C; sea-level or colder-than-absolute-zero inputs report zero, and the result is negative.
pub fn calc_dip(obs_alt_m: f64, lapse_rate_k_m: f64, atpress_mbar: f64, attemp_c: f64) -> f64 {
    if obs_alt_m <= 0.0 {
        return 0.0;
    }
    let ratio = DIP_EARTH_RADIUS_M / (DIP_EARTH_RADIUS_M + obs_alt_m);
    if ratio >= 1.0 {
        return 0.0;
    }
    let dip_geometric = ratio.acos().to_degrees();
    let krefr = (0.0342 + lapse_rate_k_m) / (0.154 * 0.0238);
    let t_kelvin = 273.16 + attemp_c;
    if t_kelvin <= 0.0 {
        return 0.0;
    }
    let mut d = 1.0 - 1.8480 * krefr * atpress_mbar / (t_kelvin * t_kelvin);
    if d < 0.0 {
        d = 0.0;
    }
    -dip_geometric * d.sqrt()
}

/// Plain `refrac` forward branch: apparent altitude for a true altitude
pub fn refrac_true_to_app(alt_deg: f64, atpress_mbar: f64, attemp_c: f64) -> f64 {
    calc_refrac_compat_true_to_app(alt_deg, atpress_mbar, attemp_c)
}

/// Plain `refrac` inverse branch: true altitude for an apparent altitude
pub fn refrac_app_to_true(alt_deg: f64, atpress_mbar: f64, attemp_c: f64) -> f64 {
    calc_refrac_compat_app_to_true(alt_deg, atpress_mbar, attemp_c)
}

/// Extended `refrac_extended` forward branch with an explicit lapse rate.
/// Returns `(converted_altitude, (true_alt, apparent_alt, refraction,
/// dip))`.
pub fn refrac_extended_true_to_app(
    alt_deg: f64,
    geoalt_m: f64,
    atpress_mbar: f64,
    attemp_c: f64,
    lapse_rate_k_m: f64,
) -> (f64, (f64, f64, f64, f64)) {
    let dip = calc_dip(geoalt_m, lapse_rate_k_m, atpress_mbar, attemp_c);
    let true_alt = alt_deg;
    let refraction =
        calc_refraction_true_to_app(alt_deg, atpress_mbar, attemp_c, geoalt_m, lapse_rate_k_m);
    let apparent_alt = true_alt + refraction;
    if apparent_alt < dip {
        return (true_alt, (true_alt, true_alt, 0.0, dip));
    }
    (apparent_alt, (true_alt, apparent_alt, refraction, dip))
}

/// Extended `refrac_extended` inverse branch with an explicit lapse rate.
/// Returns `(converted_altitude, (true_alt, apparent_alt, refraction,
/// dip))`.
pub fn refrac_extended_app_to_true(
    alt_deg: f64,
    geoalt_m: f64,
    atpress_mbar: f64,
    attemp_c: f64,
    lapse_rate_k_m: f64,
) -> (f64, (f64, f64, f64, f64)) {
    let dip = calc_dip(geoalt_m, lapse_rate_k_m, atpress_mbar, attemp_c);
    let apparent_alt = alt_deg;
    if apparent_alt < dip {
        return (apparent_alt, (apparent_alt, apparent_alt, 0.0, dip));
    }
    let refraction =
        calc_refraction_app_to_true(alt_deg, atpress_mbar, attemp_c, geoalt_m, lapse_rate_k_m);
    let true_alt = apparent_alt - refraction;
    (true_alt, (true_alt, apparent_alt, refraction, dip))
}

/// Standard-atmosphere pressure estimate from observer altitude used by
/// the `azalt` horizon entry when the caller passes `pressure <= 0`
pub fn estimate_pressure_mbar(obs_alt_m: f64) -> f64 {
    let base = 1.0 - 0.0065 * obs_alt_m / 288.0;
    STANDARD_PRESSURE_MBAR * base.max(0.0).powf(5.255)
}

// Sinclair rise/set refraction

/// Sinclair branch join in apparent-altitude degrees.
pub const SINCLAIR_BRANCH_ALT_DEG: f64 = 17.904104638432;

/// Atmospheric refraction at an apparent altitude in degrees, in degrees
///
/// Above [`SINCLAIR_BRANCH_ALT_DEG`] the `0.97 / tan(h)` high-altitude
/// asymptote (arcmin) applies; at or below it Sinclair's rational fit.
/// Bennett's pressure/temperature compensation refers to 1010 mbar /
/// 10 C; the result converts arcmin to degrees. No input guards: the
/// caller (the fixed-point core) owns finiteness, exactly like
/// the release.
/// Fit: (34.46 + 4.23*h + 0.004*h^2)/(1 + 0.505*h + 0.0845*h^2) arcmin; compensation (atpress-80)/930/(1 + 8e-5*(r+39)*(attemp-10)), result /60 to degrees.
pub fn sinclair_refraction_deg(app_alt_deg: f64, atpress_mbar: f64, attemp_c: f64) -> f64 {
    let h = app_alt_deg;
    let r = if h > SINCLAIR_BRANCH_ALT_DEG {
        0.97 / h.to_radians().tan()
    } else {
        (34.46 + 4.23 * h + 0.004 * h * h) / (1.0 + 0.505 * h + 0.0845 * h * h)
    };
    let compensated =
        (atpress_mbar - 80.0) / 930.0 / (1.0 + 0.00008 * (r + 39.0) * (attemp_c - 10.0)) * r;
    compensated / 60.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hydrostatic_exponent_matches_release_association() {
        let expected = 9.80665f64 * 0.0289644f64 / 8.314462618f64;
        assert_eq!(HYDROSTATIC_EXPONENT_K_M.to_bits(), expected.to_bits());
    }

    #[test]
    fn icao_surface_values_are_exact() {
        assert_eq!(icao_temperature_k(0.0).to_bits(), 288.15f64.to_bits());
        assert_eq!(icao_pressure_mbar(0.0).to_bits(), 1013.25f64.to_bits());
        assert_eq!(
            refractive_index(1013.25, 288.15).to_bits(),
            (1.0 + 8.060e-5 * 1013.25 / 288.15f64).to_bits()
        );
    }

    #[test]
    fn refractive_index_guards_match_release() {
        assert_eq!(refractive_index(1013.25, 0.0), 1.0);
        assert_eq!(refractive_index(1013.25, -5.0), 1.0);
        assert_eq!(refractive_index(0.0, 288.15), 1.0);
        assert_eq!(refractive_index(-3.0, 288.15), 1.0);
    }

    #[test]
    fn trace_ray_guards_match_release() {
        assert_eq!(trace_ray(0.0, 0.0, 1013.25, 15.0, 0.0065), 0.0);
        assert_eq!(trace_ray(-2.0, 0.0, 1013.25, 15.0, 0.0065), 0.0);
        assert_eq!(trace_ray(45.0, 0.0, 0.0, 15.0, 0.0065), 0.0);
        assert_eq!(trace_ray(45.0, 0.0, -4.0, 15.0, 0.0065), 0.0);
        assert_eq!(trace_ray(45.0, 0.0, 1013.25, -273.15, 0.0065), 0.0);
        assert_eq!(trace_ray(45.0, 0.0, 1013.25, -300.0, 0.0065), 0.0);
        assert_eq!(
            calc_refraction_true_to_app(30.0, 0.0, 15.0, 0.0, 0.0065),
            0.0
        );
    }

    #[test]
    fn trace_ray_extrapolation_is_linear_past_91() {
        let r90 = trace_ray(90.0, 0.0, 1013.25, 15.0, 0.0065);
        let r91 = trace_ray(91.0, 0.0, 1013.25, 15.0, 0.0065);
        let r93 = trace_ray(93.0, 0.0, 1013.25, 15.0, 0.0065);
        let slope = r91 - r90;
        let expected = r91 + slope * 2.0;
        assert!((r93 - expected).abs() < 1e-15);
        assert!(r93 >= 0.0);
    }

    #[test]
    fn compat_forward_guards_match_release() {
        assert_eq!(calc_refrac_compat_true_to_app(-5.0, 1013.25, 15.0), -5.0);
        assert_eq!(calc_refrac_compat_true_to_app(-6.0, 1013.25, 15.0), -6.0);
        assert_eq!(
            calc_refrac_compat_true_to_app(f64::INFINITY, 1013.25, 15.0),
            f64::INFINITY
        );
        let nan_out = calc_refrac_compat_true_to_app(f64::NAN, 1013.25, 15.0);
        assert!(nan_out.is_nan());
    }

    #[test]
    fn compat_forward_branches_lift_the_horizon() {
        let low = calc_refrac_compat_true_to_app(0.0, 1013.25, 15.0);
        assert!(low > 0.4 && low < 0.6);
        let high = calc_refrac_compat_true_to_app(45.0, 1013.25, 15.0);
        assert!(high > 45.0 && high < 45.05);
        let zenith = calc_refrac_compat_true_to_app(90.0, 1013.25, 15.0);
        assert!((zenith - 90.0).abs() < 1e-9);
    }

    #[test]
    fn compat_scale_zero_denominator_matches_release() {
        assert!(compat_atmosphere_scale(1013.25, -273.0).is_infinite());
        assert!(compat_atmosphere_scale(0.0, -273.0).is_nan());
        assert!(compat_atmosphere_scale(f64::NAN, -273.0).is_nan());
        let finite = compat_atmosphere_scale(1013.25, 15.0);
        assert!((finite - (1013.25 / 1010.0 * 283.0 / 288.0)).abs() < 1e-15);
    }

    #[test]
    fn dip_guards_match_release() {
        assert_eq!(calc_dip(0.0, 0.0065, 1013.25, 10.0), 0.0);
        assert_eq!(calc_dip(-50.0, 0.0065, 1013.25, 10.0), 0.0);
        assert_eq!(calc_dip(1000.0, 0.0065, 1013.25, -273.16), 0.0);
        assert_eq!(calc_dip(1000.0, 0.0065, 1013.25, -300.0), 0.0);
    }

    #[test]
    fn dip_at_1000m_is_near_0p88_deg() {
        let dip = calc_dip(1000.0, 0.0065, 1013.25, 10.0);
        assert!(dip < -0.8 && dip > -1.0);
    }

    #[test]
    fn refrac_extended_below_dipped_horizon_returns_input() {
        let (converted, details) = refrac_extended_true_to_app(-5.0, 1000.0, 1013.25, 10.0, 0.0065);
        assert_eq!(converted, -5.0);
        assert_eq!(details, (-5.0, -5.0, 0.0, details.3));
        assert!(details.3 < 0.0);
    }

    #[test]
    fn refrac_alias_and_pressure_estimate_spellings() {
        assert_eq!(
            calc_refraction_ref_true_to_app(20.0, 1013.25, 15.0, 0.0, 0.0065).to_bits(),
            calc_refraction_true_to_app(20.0, 1013.25, 15.0, 0.0, 0.0065).to_bits()
        );
        assert_eq!(
            refrac_true_to_app(10.0, 1013.25, 15.0).to_bits(),
            calc_refrac_compat_true_to_app(10.0, 1013.25, 15.0).to_bits()
        );
        let p0 = estimate_pressure_mbar(0.0);
        assert!((p0 - 1013.25).abs() < 1e-12);
        assert_eq!(estimate_pressure_mbar(100_000.0), 0.0);
    }

    #[test]
    fn profile_tables_have_release_length() {
        assert_eq!(QUAD_N, 120);
        assert_eq!(QUAD_NODE_BITS.len(), 120);
        assert_eq!(QUAD_WEIGHT_BITS.len(), 120);
        let profile = AtmosphereProfile::new(0.0, 1013.25, 288.15, 0.0065);
        assert!(profile.n_obs() > 1.0 && profile.n_obs() < 1.001);
        assert_eq!(profile.r_obs().to_bits(), 6_371_000.0f64.to_bits());
    }

    #[test]
    fn isothermal_lapse_arms_stay_finite() {
        let p = modified_pressure_mbar(5000.0, 0.0, 1013.25, 288.15, 0.0);
        assert!(p.is_finite() && p > 0.0 && p < 1013.25);
        let r = calc_refraction_true_to_app(10.0, 1013.25, 15.0, 0.0, 0.0);
        assert!(r.is_finite() && r >= 0.0);
    }

    #[test]
    fn compat_reverse_zenith_limit_matches_release_association() {
        let expected = 90.0f64 - 1e-10f64;
        assert_eq!(
            COMPAT_REVERSE_ZENITH_LIMIT_DEG.to_bits(),
            expected.to_bits()
        );
    }

    #[test]
    fn compat_inverse_guards_match_release() {
        // Non-finite inputs pass through unchanged.
        assert!(calc_refrac_compat_app_to_true(f64::NAN, 1013.25, 15.0).is_nan());
        assert_eq!(
            calc_refrac_compat_app_to_true(f64::INFINITY, 1013.25, 15.0),
            f64::INFINITY
        );
        assert_eq!(
            calc_refrac_compat_app_to_true(f64::NEG_INFINITY, 1013.25, 15.0),
            f64::NEG_INFINITY
        );
        // The exact -4.4 pole: a zero divisor returns the input (the
        // release reports this arm via ZeroDivisionError).
        assert_eq!(calc_refrac_compat_app_to_true(-4.4, 1013.25, 15.0), -4.4);
        // Zenith-limit arm: angle = app + 7.31/(app+4.4) reaches
        // 90 - 1e-10, so the input returns unchanged.
        assert_eq!(calc_refrac_compat_app_to_true(90.0, 1013.25, 15.0), 90.0);
        assert_eq!(calc_refrac_compat_app_to_true(120.0, 1013.25, 15.0), 120.0);
        // Near-pole neighbors: the release contract admits only two
        // outcomes (input unchanged, or a strictly positive true
        // altitude); both arms are pinned bitwise by the sweep.
        for app in [-4.400001, -4.399999, -4.5, -4.3] {
            let out = calc_refrac_compat_app_to_true(app, 1013.25, 15.0);
            assert!(
                out == app || (out > 0.0 && out.is_finite()),
                "near-pole {app}: unexpected {out}"
            );
        }
    }

    #[test]
    fn compat_inverse_nominal_spelling_removes_refraction() {
        // Operation-order pin: recompute with the release association.
        let app = 10.0f64;
        let scale = compat_atmosphere_scale(1013.25, 15.0);
        let angle = app + 7.31 / (app + 4.4);
        let expected = app - (1.0 / angle.to_radians().tan() * scale / 60.0);
        assert_eq!(
            calc_refrac_compat_app_to_true(app, 1013.25, 15.0).to_bits(),
            expected.to_bits()
        );
        // Physical sanity (not a numerical comparison evidence): refraction is removed.
        let out = calc_refrac_compat_app_to_true(app, 1013.25, 15.0);
        assert!(out > 0.0 && out < app);
        // Candidates that would not stay strictly above zero return the
        // input unchanged (deep-horizon arm).
        assert_eq!(calc_refrac_compat_app_to_true(-3.0, 1013.25, 15.0), -3.0);
    }

    #[test]
    fn ray_inverse_guards_match_release() {
        assert_eq!(
            calc_refraction_app_to_true(30.0, 0.0, 15.0, 0.0, 0.0065),
            0.0
        );
        assert_eq!(
            calc_refraction_app_to_true(30.0, -4.0, 15.0, 0.0, 0.0065),
            0.0
        );
        assert_eq!(
            calc_refraction_ref_app_to_true(30.0, 0.0, 15.0, 0.0, 0.0065),
            0.0
        );
    }

    #[test]
    fn ray_inverse_newton_roundtrip_sanity() {
        // Mathematical sanity of the Newton solve only (not source
        // parity — that comes from the release sweep): inverting the
        // forward-mapped altitude recovers the true altitude to the
        // iteration stop band.
        for true_alt in [0.5, 1.0, 5.0, 10.0, 30.0, 60.0, 80.0] {
            let fwd = calc_refraction_true_to_app(true_alt, 1013.25, 15.0, 0.0, 0.0065);
            let app = true_alt + fwd;
            let inv = calc_refraction_app_to_true(app, 1013.25, 15.0, 0.0, 0.0065);
            assert!(inv >= 0.0);
            assert!(
                (app - inv - true_alt).abs() < 1e-9,
                "true {true_alt}: round-trip missed by {}",
                (app - inv - true_alt).abs()
            );
        }
        // Below-horizon apparent altitudes exercise the >91 deg
        // extrapolation arm inside the Newton loop.
        let below = calc_refraction_app_to_true(-2.0, 1013.25, 15.0, 0.0, 0.0065);
        assert!(below.is_finite() && below >= 0.0);
    }

    #[test]
    fn inverse_alias_and_branch_spellings() {
        assert_eq!(
            calc_refraction_ref_app_to_true(20.0, 1013.25, 15.0, 0.0, 0.0065).to_bits(),
            calc_refraction_app_to_true(20.0, 1013.25, 15.0, 0.0, 0.0065).to_bits()
        );
        assert_eq!(
            refrac_app_to_true(10.0, 1013.25, 15.0).to_bits(),
            calc_refrac_compat_app_to_true(10.0, 1013.25, 15.0).to_bits()
        );
    }

    #[test]
    fn refrac_extended_inverse_arms_match_release() {
        // Below the dipped horizon the input returns untouched with zero
        // refraction ( `apparent_alt < dip` arm).
        let dip = calc_dip(1000.0, 0.0065, 1013.25, 10.0);
        assert!(dip < 0.0);
        let (converted, details) = refrac_extended_app_to_true(-5.0, 1000.0, 1013.25, 10.0, 0.0065);
        assert_eq!(converted, -5.0);
        assert_eq!(details, (-5.0, -5.0, 0.0, dip));
        // Nominal arm: composition identity with the kernel words.
        let (converted, details) = refrac_extended_app_to_true(10.0, 0.0, 1013.25, 15.0, 0.0065);
        let refr = calc_refraction_app_to_true(10.0, 1013.25, 15.0, 0.0, 0.0065);
        let sea_dip = calc_dip(0.0, 0.0065, 1013.25, 15.0);
        assert_eq!(converted.to_bits(), (10.0 - refr).to_bits());
        assert_eq!(details.0.to_bits(), converted.to_bits());
        assert_eq!(details.1, 10.0);
        assert_eq!(details.2.to_bits(), refr.to_bits());
        assert_eq!(details.3.to_bits(), sea_dip.to_bits());
        assert!(refr >= 0.0);
    }

    #[test]
    fn sinclair_branches_match_release_spelling() {
        // High-altitude asymptote arm ( association).
        let high = sinclair_refraction_deg(30.0, 1013.25, 15.0);
        let r = 0.97 / 30.0f64.to_radians().tan();
        let expected =
            (1013.25 - 80.0) / 930.0 / (1.0 + 0.00008 * (r + 39.0) * (15.0 - 10.0)) * r / 60.0;
        assert_eq!(high.to_bits(), expected.to_bits());
        assert!(high > 0.0 && high < 0.05);
        // Rational-fit arm at and below the join.
        let join = sinclair_refraction_deg(SINCLAIR_BRANCH_ALT_DEG, 1013.25, 15.0);
        let h = SINCLAIR_BRANCH_ALT_DEG;
        let r_fit = (34.46 + 4.23 * h + 0.004 * h * h) / (1.0 + 0.505 * h + 0.0845 * h * h);
        let expected_fit =
            (1013.25 - 80.0) / 930.0 / (1.0 + 0.00008 * (r_fit + 39.0) * (15.0 - 10.0)) * r_fit
                / 60.0;
        assert_eq!(join.to_bits(), expected_fit.to_bits());
        // Horizon refraction is ~34 arcmin at standard conditions.
        let horizon = sinclair_refraction_deg(0.0, 1013.25, 15.0);
        assert!(horizon > 0.5 && horizon < 0.65);
        // Zero apparent altitude with zero pressure still evaluates
        // (no guards on this path, like the release).
        assert!(sinclair_refraction_deg(0.0, 0.0, 15.0).is_finite());
    }
}
