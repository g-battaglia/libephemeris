//! Published mean-element evaluation and Cartesian two-body state reduction.
//!
//! These reducers do not acquire provider state or select output frames. Their
//! callers own the epoch, central mass and input-state validity contracts.

use super::{J2000_JD, MeanCoeffs};
use crate::constants;
use crate::positions::PositionResult;

/// Ascending-power polynomial coefficients `(a)`, `(e)`, `(i_deg)`,
/// `(Omega_deg | None for Earth)`, `(pi_deg)` per planet with mean
/// elements.
pub(crate) fn mean_element_coeffs(body_id: i32) -> Option<MeanCoeffs> {
    match body_id {
        constants::MERCURY => Some(MeanCoeffs {
            a: &[0.387098310],
            e: &[0.20563175, 0.000020407, -0.0000000283, -0.00000000018],
            i: &[7.004986, 0.0018215, -0.00001810, 0.000000056],
            omega: Some(&[48.330893, 1.1861890, 0.00017587, 0.000000211]),
            pi: &[77.456119, 1.5564775, 0.00029589, 0.000000056],
        }),
        constants::VENUS => Some(MeanCoeffs {
            a: &[0.723329820],
            e: &[0.00677188, -0.000047766, 0.0000000975, 0.00000000044],
            i: &[3.394662, 0.0010037, -0.00000088, -0.000000007],
            omega: Some(&[76.679920, 0.9011190, 0.00040665, -0.000000080]),
            pi: &[131.563707, 1.4022188, -0.00107337, -0.000005315],
        }),
        constants::EARTH => Some(MeanCoeffs {
            a: &[1.000001018],
            e: &[0.01670862, -0.000042037, -0.0000001236, 0.00000000004],
            i: &[0.0],
            omega: None,
            pi: &[102.937348, 1.7195269, 0.00045962, 0.000000499],
        }),
        constants::MARS => Some(MeanCoeffs {
            a: &[1.523679342],
            e: &[0.09340062, 0.000090483, -0.0000000806, -0.00000000035],
            i: &[1.849726, -0.0006010, 0.00001276, -0.000000006],
            omega: Some(&[49.558093, 0.7720923, 0.00001605, 0.000002325]),
            pi: &[336.060234, 1.8410331, 0.00013515, 0.000000318],
        }),
        constants::JUPITER => Some(MeanCoeffs {
            a: &[5.202603191, 0.0000001913],
            e: &[0.04849485, 0.000163244, -0.0000004719, -0.00000000197],
            i: &[1.303270, -0.0054966, 0.00000465, -0.000000004],
            omega: Some(&[100.464441, 1.0209550, 0.00040117, 0.000000569]),
            pi: &[14.331309, 1.6126668, 0.00103127, -0.000004569],
        }),
        constants::SATURN => Some(MeanCoeffs {
            a: &[9.554909596, -0.0000021389],
            e: &[0.05550862, -0.000346818, -0.0000006456, 0.00000000338],
            i: &[2.488878, -0.0037363, -0.00001516, 0.000000089],
            omega: Some(&[113.665524, 0.8770979, -0.00012067, -0.000002380]),
            pi: &[93.056787, 1.9637694, 0.00083757, 0.000004899],
        }),
        constants::URANUS => Some(MeanCoeffs {
            a: &[19.218446062, -0.0000000372, 0.00000000098],
            e: &[0.04629590, -0.000027337, 0.0000000790, 0.00000000025],
            i: &[0.773196, 0.0007744, 0.00003749, -0.000000092],
            omega: Some(&[74.005947, 0.5211258, 0.00133982, 0.000018516]),
            pi: &[173.005159, 1.4863784, 0.00021450, 0.000000433],
        }),
        constants::NEPTUNE => Some(MeanCoeffs {
            a: &[30.110386869, -0.0000001663, 0.00000000069],
            e: &[0.00898809, 0.000006408, -0.0000000008, -0.00000000005],
            i: &[1.769952, -0.0093082, -0.00000708, 0.000000028],
            omega: Some(&[131.784057, 1.1022057, 0.00026006, -0.000000636]),
            pi: &[48.123691, 1.4262677, 0.00037918, -0.000000003],
        }),
        _ => None,
    }
}

/// Ascending-power polynomial evaluation.
pub(crate) fn poly_asc(coeffs: &[f64], t: f64) -> f64 {
    let mut result = 0.0;
    for c in coeffs.iter().rev() {
        result = result * t + c;
    }
    result
}

/// Mean orbital elements at the mean equinox of date: `(a_au, e, i_deg, Omega_deg, pi_deg)`.
/// Earth reports `Omega = 0` and `i = 0` (no node); bodies without mean
/// elements (Pluto and everything else) return `None` so callers fall
/// back to osculating fits. Non-finite epochs return `None` (fail-closed;
/// the model would propagate NaN through the polynomial). Coefficients
/// transcribe Meeus, Astronomical Algorithms (2 nd ed.), Table 31.A, from
/// Simon et al. (1994, A&A 282, 663-683); `T` is Julian centuries of 36525
/// days TT from J2000, referred to the mean ecliptic and equinox of date.
/// `Omega` and `pi` are normalized to [0, 360) degrees.
pub fn mean_orbital_elements(body_id: i32, jd_tt: f64) -> Option<[f64; 5]> {
    let el = mean_element_coeffs(body_id)?;
    if !jd_tt.is_finite() {
        return None;
    }
    let t = (jd_tt - J2000_JD) / 36525.0;
    let omega_deg = match el.omega {
        Some(coeffs) => poly_asc(coeffs, t).rem_euclid(360.0),
        None => 0.0,
    };
    Some([
        poly_asc(el.a, t),
        poly_asc(el.e, t),
        poly_asc(el.i, t),
        omega_deg,
        poly_asc(el.pi, t).rem_euclid(360.0),
    ])
}

/// Two-body semi-major axis from a Cartesian state and gravitational parameter.
///
/// Position components are AU, velocity components AU/day and `gm` AU³/day²
/// about the same center. Returns AU, negative for a hyperbolic orbit. A zero
/// or nonfinite position norm, or an exactly zero energy denominator, returns
/// `0.0`. Other invalid inputs or overflowing arithmetic can produce NaN or
/// infinity; this unchecked reducer is not a general state-domain validator.
pub fn osculating_semi_major_axis(
    rx: f64,
    ry: f64,
    rz: f64,
    vx: f64,
    vy: f64,
    vz: f64,
    gm: f64,
) -> f64 {
    let r = (rx * rx + ry * ry + rz * rz).sqrt();
    if !r.is_finite() || r <= 0.0 {
        return 0.0;
    }
    let v2 = vx * vx + vy * vy + vz * vz;
    let denom = 2.0 / r - v2 / gm;
    if denom != 0.0 { 1.0 / denom } else { 0.0 }
}

/// Ecliptic spherical state to Cartesian state in the same frame
/// : `(lon, lat, r,
/// dlon, dlat, ddist)` in degrees/degrees-per-day/AU to `(x, y, z, vx,
/// vy, vz)` in AU/AU-day, with the analytic time derivative so the pair
/// is self-consistent for a two-body reduction.
pub fn ecliptic_spherical_to_cartesian_state(pos: PositionResult) -> [f64; 6] {
    let [lon, lat, r, dlon, dlat, ddist] = pos;
    let lam = lon.to_radians();
    let bet = lat.to_radians();
    let dlam = dlon.to_radians();
    let dbet = dlat.to_radians();
    let cos_lam = lam.cos();
    let sin_lam = lam.sin();
    let cos_bet = bet.cos();
    let sin_bet = bet.sin();
    [
        r * cos_bet * cos_lam,
        r * cos_bet * sin_lam,
        r * sin_bet,
        ddist * cos_bet * cos_lam - r * sin_bet * cos_lam * dbet - r * cos_bet * sin_lam * dlam,
        ddist * cos_bet * sin_lam - r * sin_bet * sin_lam * dbet + r * cos_bet * cos_lam * dlam,
        ddist * sin_bet + r * cos_bet * dbet,
    ]
}

/// Reduce an ecliptic state vector to the 50-slot osculating-element
/// tuple.
///
/// `state` is `(x, y, z, vx, vy, vz)` in AU and AU/day. `gm` is the
/// gravitational parameter driving the two-body geometry, in AU³/day²;
/// `gm_slots` optionally selects a different central-mass parameter in the
/// same units for period/rate slots 10–14 (`None` keeps `gm` everywhere).
/// Angular momentum is `h = r × v`, orbital energy is `|v|²/2 - gm/|r|`, and
/// the eccentricity vector is `((|v|²/gm - 1/|r|)r - (r·v/gm)v)`.
/// `tt_jd` is a finite JD TT for the perihelion epoch; `body_id` selects the
/// Earth's zero synodic-period special case and does not select a mass model.
///
/// For positive finite `a`, the period parameter `μ` is `gm_slots` when
/// supplied, otherwise `gm`, and `n = sqrt(μ/a³)` is radians/day. The sidereal
/// period is `P = 2π/n` days, reported as `P/365.24219` years. With `n_deg`
/// in degrees/day and precession `p = 50.2882/3600/365.25` degrees/day, the
/// tropical period is `360/(n_deg + p)/365.24219` years. Perihelion passage is
/// `T = tt_jd - M/n`, with mean anomaly `M` in radians. Non-positive or
/// non-finite `a` sets `n = 0`; absent positive motion yields infinite period
/// years and a zero perihelion epoch.
///
/// The signed synodic period uses an Earth sidereal year `P_E = 365.256363`
/// days: `P_syn = P_E*P/(P - P_E)`. It is zero for Earth or an absent positive
/// finite `P`, and infinite when `abs(P - P_E) <= 1e-10` days. Interior
/// orbits can therefore retain a negative synodic period.
///
/// | Slot | Element | Unit |
/// | --- | --- | --- |
/// | 0 | Semi-major axis `a` | AU |
/// | 1 | Eccentricity `e` | Dimensionless |
/// | 2 | Inclination `i` | Degrees |
/// | 3 | Ascending-node longitude `Ω` | Degrees |
/// | 4 | Perihelion argument `ω` | Degrees |
/// | 5 | Perihelion longitude `ϖ` | Degrees |
/// | 6 | Mean anomaly `M` | Degrees |
/// | 7 | True anomaly `ν` | Degrees |
/// | 8 | Eccentric anomaly `E` | Degrees |
/// | 9 | Mean longitude `L` | Degrees |
/// | 10 | Sidereal period | Years of 365.24219 days |
/// | 11 | Mean daily motion `n` | Degrees/day |
/// | 12 | Tropical period | Years of 365.24219 days |
/// | 13 | Signed synodic period | Days |
/// | 14 | Perihelion-passage epoch | JD TT |
/// | 15 | Perihelion distance `q` | AU |
/// | 16 | Aphelion distance `Q` | AU |
/// | 17–49 | Reserved | Positive zero |
///
/// Scalar norms and energy retain floating-point powers, rather than replacing
/// them with repeated multiplication: rounding near an apsis can decide the
/// anomaly's acos argument. Operation order
/// follows the model, including the near-circular
/// (`e < 1e-10`) and near-equatorial (`n_mag <= 1e-10`) degenerate arms
/// and the parabolic/hyperbolic guards. The `% (2*pi)` folds act on
/// non-negative sums, where truncated and Euclidean remainders
/// (floored) agree. Degenerate guards cover near-parabolic energy
/// (`|eps|` <= 1 e-10), near-circular (`e` < 1 e-10) and near-equatorial
/// (`n_mag` <= 1 e-10) arms. This unchecked numerical reducer requires finite
/// state components, finite `tt_jd`, a positive finite position norm and
/// positive finite `gm`
/// (and `gm_slots` when supplied); intermediate norms/powers must not overflow.
/// Callers are responsible for rejecting unusable states before this reduction.
/// Invalid inputs are not converted into a typed error here; arithmetic may
/// produce NaN or infinity outside this reducer's supported domain.
pub fn orbital_elements_from_ecliptic_state(
    state: [f64; 6],
    gm: f64,
    body_id: i32,
    tt_jd: f64,
    gm_slots: Option<f64>,
) -> [f64; 50] {
    let tau = 2.0 * std::f64::consts::PI;
    let [x, y, z, vx, vy, vz] = state;
    let r_mag = (x.powf(2.0) + y.powf(2.0) + z.powf(2.0)).sqrt();
    let v_mag = (vx.powf(2.0) + vy.powf(2.0) + vz.powf(2.0)).sqrt();
    let r_dot_v = x * vx + y * vy + z * vz;
    let hx = y * vz - z * vy;
    let hy = z * vx - x * vz;
    let hz = x * vy - y * vx;
    let h_mag = (hx.powf(2.0) + hy.powf(2.0) + hz.powf(2.0)).sqrt();
    let nx = -hy;
    let ny = hx;
    let n_mag = (nx.powf(2.0) + ny.powf(2.0)).sqrt();
    let e_coef = v_mag.powf(2.0) / gm - 1.0 / r_mag;
    let rv_coef = r_dot_v / gm;
    let ex = e_coef * x - rv_coef * vx;
    let ey = e_coef * y - rv_coef * vy;
    let ez = e_coef * z - rv_coef * vz;
    let e_mag = (ex.powf(2.0) + ey.powf(2.0) + ez.powf(2.0)).sqrt();
    let epsilon = v_mag.powf(2.0) / 2.0 - gm / r_mag;
    let a = if epsilon.abs() > 1e-10 {
        -gm / (2.0 * epsilon)
    } else {
        f64::INFINITY
    };
    let mut e = e_mag;
    if e < 1e-10 {
        e = 0.0;
    }
    let incl = if h_mag > 1e-10 {
        (hz / h_mag).clamp(-1.0, 1.0).acos()
    } else {
        0.0
    };
    let omega_node = if n_mag > 1e-10 {
        let mut o = ny.atan2(nx);
        if o < 0.0 {
            o += tau;
        }
        o
    } else {
        0.0
    };
    let omega_peri = if n_mag > 1e-10 && e_mag > 1e-10 {
        let n_dot_e = nx * ex + ny * ey;
        let cos_omega = (n_dot_e / (n_mag * e_mag)).clamp(-1.0, 1.0);
        let mut o = cos_omega.acos();
        if ez < 0.0 {
            o = tau - o;
        }
        o
    } else if e_mag > 1e-10 {
        let mut o = ey.atan2(ex);
        if o < 0.0 {
            o += tau;
        }
        o
    } else {
        0.0
    };
    let nu = if e_mag > 1e-10 {
        let e_dot_r = ex * x + ey * y + ez * z;
        let cos_nu = (e_dot_r / (e_mag * r_mag)).clamp(-1.0, 1.0);
        let mut n = cos_nu.acos();
        if r_dot_v < 0.0 {
            n = tau - n;
        }
        n
    } else if n_mag > 1e-10 {
        let n_dot_r = nx * x + ny * y;
        let cos_nu = (n_dot_r / (n_mag * r_mag)).clamp(-1.0, 1.0);
        let mut n = cos_nu.acos();
        if z < 0.0 {
            n = tau - n;
        }
        n
    } else {
        let mut n = y.atan2(x);
        if n < 0.0 {
            n += tau;
        }
        n
    };
    let ecc_anom = if e < 1.0 && e > 1e-10 {
        let tan_half_nu = (nu / 2.0).tan();
        let sqrt_term = ((1.0 - e) / (1.0 + e)).sqrt();
        let mut ea = 2.0 * (sqrt_term * tan_half_nu).atan();
        if ea < 0.0 {
            ea += tau;
        }
        ea
    } else {
        nu
    };
    let mean_anom = if e < 1.0 {
        let mut m = ecc_anom - e * ecc_anom.sin();
        if m < 0.0 {
            m += tau;
        }
        m
    } else {
        ecc_anom
    };
    let mean_long = (omega_node + omega_peri + mean_anom) % tau;
    let varpi = (omega_node + omega_peri) % tau;
    // The period/rate slots use the central mass alone (compatibility
    // contract); the geometry above already used the physical two-body
    // `gm`. Retain the scalar power's rounding in the period calculation.
    let gm_for_slots = gm_slots.unwrap_or(gm);
    let n = if a > 0.0 && a < f64::INFINITY {
        (gm_for_slots / a.powf(3.0)).sqrt()
    } else {
        0.0
    };
    let (q, qq) = if a > 0.0 && a < f64::INFINITY {
        (a * (1.0 - e), a * (1.0 + e))
    } else {
        (r_mag, f64::INFINITY)
    };
    let (p_days, p_sid_years) = if n > 0.0 {
        let days = tau / n;
        (days, days / 365.24219)
    } else {
        (f64::INFINITY, f64::INFINITY)
    };
    let n_deg = n.to_degrees();
    // Tropical period via the IAU general precession in longitude
    // (~50.29"/yr); the external sidereal/tropical-year factor some
    // implementations carry has no published basis and is intentionally
    // not reproduced.
    let precession_deg_per_day = 50.2882 / 3600.0 / 365.25;
    let p_trop_years = if n_deg > 0.0 {
        360.0 / (n_deg + precession_deg_per_day) / 365.24219
    } else {
        f64::INFINITY
    };
    let t_jd = if n > 0.0 { tt_jd - mean_anom / n } else { 0.0 };
    // Synodic period from 1/P_syn = 1/P_E - 1/P between sidereal
    // periods; Earth enters as the sidereal year (365.256363 d).
    let p_earth_days = 365.256363;
    let p_syn = if p_days > 0.0 && p_days < f64::INFINITY && body_id != constants::EARTH {
        let denom = p_days - p_earth_days;
        if denom.abs() > 1e-10 {
            p_earth_days * p_days / denom
        } else {
            f64::INFINITY
        }
    } else {
        0.0
    };
    let mut out = [0.0; 50];
    out[0] = a;
    out[1] = e;
    out[2] = incl.to_degrees();
    out[3] = omega_node.to_degrees();
    out[4] = omega_peri.to_degrees();
    out[5] = varpi.to_degrees();
    out[6] = mean_anom.to_degrees();
    out[7] = nu.to_degrees();
    out[8] = ecc_anom.to_degrees();
    out[9] = mean_long.to_degrees();
    out[10] = p_sid_years;
    out[11] = n_deg;
    out[12] = p_trop_years;
    out[13] = p_syn;
    out[14] = t_jd;
    out[15] = q;
    out[16] = qq;
    out
}
