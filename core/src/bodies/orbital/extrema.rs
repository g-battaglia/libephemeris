//! Distance extrema of fixed osculating ellipses.
//!
//! The coarse scan and local windows deliberately retain their iteration and
//! strict-comparison order: equal distances keep the first sampled pair.

use super::{
    ORBIT_EXTREMUM_COARSE_N, ORBIT_EXTREMUM_REFINE_ITERS, ORBIT_EXTREMUM_REFINE_WIDTH, OrbitEllipse,
};
use crate::errors::{Error, Result};

impl OrbitEllipse {
    /// Store finite elements with a non-negative semi-major axis in AU.
    ///
    /// Eccentricity is dimensionless; inclination, node and perihelion argument
    /// are radians in one common reference frame. Nonfinite fields or a negative
    /// axis return `NonFiniteArgument`. Zero axes and eccentricities outside
    /// `[0, 1)` are admitted: this constructor does not certify an ellipse.
    /// Physical elliptic searches require the caller to choose usable elements.
    pub fn new(a_au: f64, e: f64, inc_rad: f64, node_rad: f64, arg_peri_rad: f64) -> Result<Self> {
        for v in [a_au, e, inc_rad, node_rad, arg_peri_rad] {
            if !v.is_finite() {
                return Err(Error::NonFiniteArgument { value: v });
            }
        }
        if a_au < 0.0 {
            return Err(Error::NonFiniteArgument { value: a_au });
        }
        Ok(Self {
            a_au,
            e,
            inc_rad,
            node_rad,
            arg_peri_rad,
        })
    }

    /// Semi-major axis in AU.
    pub fn a_au(self) -> f64 {
        self.a_au
    }

    /// Eccentricity.
    pub fn e(self) -> f64 {
        self.e
    }

    /// Inclination in radians.
    pub fn inc_rad(self) -> f64 {
        self.inc_rad
    }

    /// Longitude of ascending node in radians.
    pub fn node_rad(self) -> f64 {
        self.node_rad
    }

    /// Argument of perihelion in radians.
    pub fn arg_peri_rad(self) -> f64 {
        self.arg_peri_rad
    }
}

/// Read the extremum-search ellipse from a 50-slot osculating-element tuple.
///
/// Slots 0–4 are `a` in AU, eccentricity, inclination, node and perihelion
/// argument in degrees, as in [`super::orbital_elements_from_ecliptic_state`].
/// Angles convert to radians; admission follows [`OrbitEllipse::new`]. The
/// input tuple, not this converter, determines the center and reference plane.
pub fn orbit_ellipse_from_elements(el: &[f64; 50]) -> Result<OrbitEllipse> {
    OrbitEllipse::new(
        el[0],
        el[1],
        el[2].to_radians(),
        el[3].to_radians(),
        el[4].to_radians(),
    )
}

/// Ecliptic Cartesian position (AU) at eccentric anomaly `anom_rad`
/// (scalar form of the vector kernel:
/// perifocal plane with the `max(0, 1 - e^2)` floor, rotation by
/// argument of perihelion, inclination about the node line, rotation by
/// the ascending node, in that order).
pub fn ellipse_position_au(el: OrbitEllipse, anom_rad: f64) -> [f64; 3] {
    let cos_e = anom_rad.cos();
    let sin_e = anom_rad.sin();
    let x_pf = el.a_au * (cos_e - el.e);
    let y_pf = el.a_au * (1.0 - el.e * el.e).max(0.0).sqrt() * sin_e;
    let cos_w = el.arg_peri_rad.cos();
    let sin_w = el.arg_peri_rad.sin();
    let x1 = x_pf * cos_w - y_pf * sin_w;
    let y1 = x_pf * sin_w + y_pf * cos_w;
    let cos_i = el.inc_rad.cos();
    let sin_i = el.inc_rad.sin();
    let y2 = y1 * cos_i;
    let z2 = y1 * sin_i;
    let cos_o = el.node_rad.cos();
    let sin_o = el.node_rad.sin();
    [x1 * cos_o - y2 * sin_o, x1 * sin_o + y2 * cos_o, z2]
}

/// Squared distance between two ellipse positions with the model
/// `maximum(d_sq, 0.0)` floor (never a NaN from rounding below zero).
fn ellipse_pair_d_sq(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).max(0.0)
}

/// Extremum of `|P1(E1) - P2(E2)|` over one anomaly window pair sampled
/// with `ORBIT_EXTREMUM_REFINE_WIDTH` points per axis.
/// Returns the winning center pair and its squared distance.
fn extremum_window(
    el1: OrbitEllipse,
    el2: OrbitEllipse,
    c1: f64,
    c2: f64,
    half: f64,
    want_max: bool,
) -> ((f64, f64), f64) {
    let n = ORBIT_EXTREMUM_REFINE_WIDTH;
    let mut best_d_sq = if want_max {
        f64::NEG_INFINITY
    } else {
        f64::INFINITY
    };
    let mut best = (c1, c2);
    for i in 0..n {
        let g1 = if n == 1 {
            c1
        } else {
            c1 - half + 2.0 * half * (i as f64) / ((n - 1) as f64)
        };
        let p1 = ellipse_position_au(el1, g1);
        for j in 0..n {
            let g2 = if n == 1 {
                c2
            } else {
                c2 - half + 2.0 * half * (j as f64) / ((n - 1) as f64)
            };
            let d_sq = ellipse_pair_d_sq(p1, ellipse_position_au(el2, g2));
            if (want_max && d_sq > best_d_sq) || (!want_max && d_sq < best_d_sq) {
                best_d_sq = d_sq;
                best = (g1, g2);
            }
        }
    }
    (best, best_d_sq)
}

/// Approximate 3D maximum and minimum distances between two osculating ellipses.
///
/// A fixed anomaly grid selects the best sampled cell for each extremum. A
/// local shrinking window (`half *= 0.25`, `ORBIT_EXTREMUM_REFINE_ITERS`
/// iterations) refines only those cells. This bounded search supplies neither
/// a global-optimality certificate nor a uniform absolute-error bound.
/// Returns `(max_distance_au, min_distance_au)`.
///
/// Each sampled pair uses component differences `dx`, `dy`, `dz` and
/// `(dx * dx + dy * dy + dz * dz).max(0.0)` for squared distance, in that
/// order. The square root applies once to each final refined extremum.
pub fn two_ellipse_min_max(el1: OrbitEllipse, el2: OrbitEllipse) -> (f64, f64) {
    let two_pi = 2.0 * std::f64::consts::PI;
    let n = ORBIT_EXTREMUM_COARSE_N;
    let mut coarse: Vec<[f64; 3]> = Vec::with_capacity(n);
    for k in 0..n {
        coarse.push(ellipse_position_au(el1, two_pi * (k as f64) / (n as f64)));
    }
    let mut coarse2: Vec<[f64; 3]> = Vec::with_capacity(n);
    for k in 0..n {
        coarse2.push(ellipse_position_au(el2, two_pi * (k as f64) / (n as f64)));
    }
    let (mut i_min, mut i_max) = ((0usize, 0usize), (0usize, 0usize));
    let (mut d_min, mut d_max) = (f64::INFINITY, f64::NEG_INFINITY);
    for (i, p1) in coarse.iter().enumerate() {
        for (j, p2) in coarse2.iter().enumerate() {
            let d_sq = ellipse_pair_d_sq(*p1, *p2);
            if d_sq < d_min {
                d_min = d_sq;
                i_min = (i, j);
            }
            if d_sq > d_max {
                d_max = d_sq;
                i_max = (i, j);
            }
        }
    }
    let refine = |c1: f64, c2: f64, want_max: bool| -> f64 {
        let (mut cc1, mut cc2) = (c1, c2);
        let mut half = two_pi / (n as f64);
        let mut best = 0.0;
        for _ in 0..ORBIT_EXTREMUM_REFINE_ITERS {
            let ((n1, n2), d_sq) = extremum_window(el1, el2, cc1, cc2, half, want_max);
            cc1 = n1;
            cc2 = n2;
            best = d_sq;
            half *= 0.25;
        }
        best.sqrt()
    };
    let grid = |k: usize| -> f64 { two_pi * (k as f64) / (n as f64) };
    let d_max_r = refine(grid(i_max.0), grid(i_max.1), true);
    let d_min_r = refine(grid(i_min.0), grid(i_min.1), false);
    (d_max_r, d_min_r)
}

/// Heliocentric/barycentric extremes from one element tuple: `None` when `a <= 0`
/// . Otherwise `(aphelion, perihelion)` (slots
/// `[16], [15]`, i.e. the `(max, min)` order).
pub fn helio_distance_extremes(el: &[f64; 50]) -> Option<(f64, f64)> {
    if el[0] <= 0.0 {
        return None;
    }
    Some((el[16], el[15]))
}

/// Geocentric planet/asteroid extremes between the body's ellipse and
/// the EMB ellipse: `None` for an
/// invalid or non-elliptical osculating orbit (`a <= 0` or `e >= 1`,
/// where the model returns the zero triple). Otherwise returns the 3D
/// `(max, min)` in AU from the bounded approximation in [`two_ellipse_min_max`].
/// The two conics must share heliocentric axes and units. Other nonfinite
/// element fields can return constructor errors; the early non-elliptic gate
/// is not a complete validator. The result is not a simultaneous encounter
/// distance, trajectory extremum or certified global optimum.
pub fn geocentric_two_ellipse_extremes(
    el_body: &[f64; 50],
    el_emb: &[f64; 50],
) -> Result<Option<(f64, f64)>> {
    if el_body[0] <= 0.0 || el_body[1] >= 1.0 {
        return Ok(None);
    }
    let b = orbit_ellipse_from_elements(el_body)?;
    let e = orbit_ellipse_from_elements(el_emb)?;
    Ok(Some(two_ellipse_min_max(b, e)))
}

/// Geocentric-Sun extremes from the EMB elements:
/// the EMB perihelion/aphelion `(max, min)`.
pub fn sun_extremes_from_emb_elements(el_emb: &[f64; 50]) -> (f64, f64) {
    (el_emb[16], el_emb[15])
}

/// Geocentric-Earth extremes from the EMB elements:
/// coincident minimum `0` and orbit-diameter maximum `2a`.
pub fn earth_extremes_from_emb_elements(el_emb: &[f64; 50]) -> (f64, f64) {
    (2.0 * el_emb[0], 0.0)
}
