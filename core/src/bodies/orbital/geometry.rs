//! Orbit-frame conics and lunar node/apse composition.
//!
//! Geometry consumes explicit elements or states; provider and observer
//! acquisition belongs to the serving layer.

use super::{OrbitShape, nodaps_dir_to_j2000};
use crate::errors::{Error, Result};
use crate::positions::PositionResult;

impl OrbitShape {
    /// Semi-major axis (AU).
    pub fn semi_major_axis(self) -> f64 {
        self.a
    }
    /// Eccentricity.
    pub fn eccentricity(self) -> f64 {
        self.e
    }
    /// Inclination (radians).
    pub fn inclination(self) -> f64 {
        self.incl
    }
    /// Longitude of the ascending node (radians).
    pub fn node(self) -> f64 {
        self.node
    }
    /// Argument of perihelion (radians).
    pub fn arg_peri(self) -> f64 {
        self.arg_peri
    }
    /// Semi-latus rectum (AU).
    pub fn semi_latus_rectum(self) -> f64 {
        self.p
    }
}

/// Orbit shape from mean elements `(a, e, i_deg, Om_deg, pi_deg)`
/// : `omega = pi - Om`, folded
/// through Euclidean modulo (`rem_euclid` for the
/// possibly-negative difference). The shape lives in the
/// mean-equinox-of-date frame of the published mean elements (Meeus
/// Table 31.A / Simon et al. 1994); aligning it with true-equinox-of-date
/// output takes the nutation-in-longitude rotation applied at the serving
/// layer, not here.
pub fn mean_shape_from_elements(
    a: f64,
    e: f64,
    i_deg: f64,
    om_deg: f64,
    pi_deg: f64,
) -> OrbitShape {
    OrbitShape {
        a,
        e,
        incl: i_deg.to_radians(),
        node: om_deg.to_radians(),
        arg_peri: (pi_deg - om_deg).rem_euclid(360.0).to_radians(),
        p: a * (1.0 - e * e),
    }
}

/// Orbit shape from an ecliptic state vector and central `gm`
/// : angular momentum,
/// node vector, eccentricity vector, argument of perihelion and the
/// semi-latus rectum, with the model `1e-10` degenerate guards.
///
/// Position is AU, velocity AU/day and `gm` AU³/day² about the same center
/// and in the same ecliptic axes at one epoch. Returned shape angles are
/// radians and lengths AU. The caller supplies finite, usable two-body states
/// and positive `gm`; this unchecked reducer can produce NaN or infinity for
/// singular states and does not return a typed domain error.
pub fn osculating_shape_from_ecliptic(r_ecl: [f64; 3], v_ecl: [f64; 3], gm: f64) -> OrbitShape {
    let r_mag = (r_ecl[0] * r_ecl[0] + r_ecl[1] * r_ecl[1] + r_ecl[2] * r_ecl[2]).sqrt();
    let v_mag = (v_ecl[0] * v_ecl[0] + v_ecl[1] * v_ecl[1] + v_ecl[2] * v_ecl[2]).sqrt();
    let hx = r_ecl[1] * v_ecl[2] - r_ecl[2] * v_ecl[1];
    let hy = r_ecl[2] * v_ecl[0] - r_ecl[0] * v_ecl[2];
    let hz = r_ecl[0] * v_ecl[1] - r_ecl[1] * v_ecl[0];
    let h_mag = (hx * hx + hy * hy + hz * hz).sqrt();
    let incl = if h_mag > 0.0 {
        (hz / h_mag).clamp(-1.0, 1.0).acos()
    } else {
        0.0
    };
    let nx = -hy;
    let ny = hx;
    let n_mag = (nx * nx + ny * ny).sqrt();
    let tau = 2.0 * std::f64::consts::PI;
    let node = if n_mag > 1e-10 {
        let mut o = ny.atan2(nx);
        if o < 0.0 {
            o += tau;
        }
        o
    } else {
        0.0
    };
    let r_dot_v = r_ecl[0] * v_ecl[0] + r_ecl[1] * v_ecl[1] + r_ecl[2] * v_ecl[2];
    let coef1 = v_mag * v_mag / gm - 1.0 / r_mag;
    let coef2 = r_dot_v / gm;
    let ex = coef1 * r_ecl[0] - coef2 * v_ecl[0];
    let ey = coef1 * r_ecl[1] - coef2 * v_ecl[1];
    let ez = coef1 * r_ecl[2] - coef2 * v_ecl[2];
    let e_mag = (ex * ex + ey * ey + ez * ez).sqrt();
    let arg_peri = if n_mag > 1e-10 && e_mag > 1e-10 {
        let n_dot_e = nx * ex + ny * ey;
        let cos_omega = (n_dot_e / (n_mag * e_mag)).clamp(-1.0, 1.0);
        let mut o = cos_omega.acos();
        if ez < 0.0 {
            o = tau - o;
        }
        o
    } else {
        0.0
    };
    let a = 1.0 / (2.0 / r_mag - v_mag * v_mag / gm);
    let p = if e_mag < 1.0 {
        a * (1.0 - e_mag * e_mag)
    } else {
        h_mag * h_mag / gm
    };
    OrbitShape {
        a,
        e: e_mag,
        incl,
        node,
        arg_peri,
        p,
    }
}

/// Ecliptic position at true anomaly `nu` in radians.
///
/// Return `(position_au, radius_au)` relative to the conic's primary focus in
/// the reference axes of `shape`. The conic radius uses the `1e-10` denominator
/// guard, then the 3-1-3 perifocal chain `(node, incl, arg_peri)` rotates it.
/// A near-zero denominator substitutes the stored semi-major axis; this is
/// not a physical point at infinity for a parabolic/hyperbolic asymptote.
/// No input validation or provider acquisition occurs here.
pub fn orbit_point_3d(shape: OrbitShape, nu: f64) -> ([f64; 3], f64) {
    let denom = 1.0 + shape.e * nu.cos();
    let r_orb = if denom.abs() > 1e-10 {
        shape.p / denom
    } else {
        shape.a
    };
    let x_orb = r_orb * nu.cos();
    let y_orb = r_orb * nu.sin();
    let (sin_o, cos_o) = shape.node.sin_cos();
    let (sin_i, cos_i) = shape.incl.sin_cos();
    let (sin_w, cos_w) = shape.arg_peri.sin_cos();
    let xe = (cos_o * cos_w - sin_o * sin_w * cos_i) * x_orb
        + (-cos_o * sin_w - sin_o * cos_w * cos_i) * y_orb;
    let ye = (sin_o * cos_w + cos_o * sin_w * cos_i) * x_orb
        + (-sin_o * sin_w + cos_o * cos_w * cos_i) * y_orb;
    let ze = (sin_w * sin_i) * x_orb + (cos_w * sin_i) * y_orb;
    ([xe, ye, ze], r_orb)
}

/// Intersect the conic with its reference ecliptic plane.
///
/// At the ascending/descending node the argument of latitude is exactly
/// zero/π. The radius is therefore `p/(1 ± e cos(ω))` and the unit direction
/// is `±(cos(Ω), sin(Ω), 0)`. Evaluating this identity directly avoids a
/// cancellation residual perpendicular to the defining plane. Units and
/// degenerate-denominator behavior match [`orbit_point_3d`].
pub(super) fn orbit_node_3d(shape: OrbitShape, ascending: bool) -> ([f64; 3], f64) {
    let direction = if ascending { 1.0 } else { -1.0 };
    let denominator = 1.0 + direction * shape.e * shape.arg_peri.cos();
    let radius = if denominator.abs() > 1e-10 {
        shape.p / denominator
    } else {
        shape.a
    };
    let (sin_node, cos_node) = shape.node.sin_cos();
    (
        [
            direction * radius * cos_node,
            direction * radius * sin_node,
            0.0,
        ],
        radius,
    )
}

/// Second (empty) focal point of the ellipse:
/// `2ae` from the primary focus along the apse line in the
/// anti-perihelion direction. Carried in the aphelion slot under
/// `NODBIT_FOPOINT`.
pub fn focal_point_3d(shape: OrbitShape) -> ([f64; 3], f64) {
    let f_dist = 2.0 * shape.a * shape.e;
    let fx_orb = -f_dist;
    let (sin_o, cos_o) = shape.node.sin_cos();
    let (sin_i, cos_i) = shape.incl.sin_cos();
    let (sin_w, cos_w) = shape.arg_peri.sin_cos();
    let fxe = (cos_o * cos_w - sin_o * sin_w * cos_i) * fx_orb;
    let fye = (sin_o * cos_w + cos_o * sin_w * cos_i) * fx_orb;
    let fze = (sin_w * sin_i) * fx_orb;
    ([fxe, fye, fze], f_dist)
}

/// Assemble the four lunar node/apse points from lunar-theory places
/// : node distances from the conic
/// radius at the node's true anomaly, perigee opposite the apogee, the
/// `NODBIT_FOPOINT` second-focus substitution, heliocentric/
/// barycentric re-referencing (helio wins when both are set), the J2000
/// direction rule, and the optional mean-ayanamsha subtraction
/// (`None` keeps tropical longitudes; the serving refuses sidereal
/// until the sidereal serving layer binds the ayanamsha). Errors: returns
/// `Error::NonFiniteArgument` on any non-finite angle, distance, or `jd_tt`
/// input, or on a non-positive apogee/perigee distance (a
/// checked-construction addition; the model branch itself is total on
/// finite inputs).
#[allow(clippy::too_many_arguments)]
pub fn moon_nodaps_points(
    node_lon: f64,
    node_lat: f64,
    apog_lon: f64,
    apog_lat: f64,
    apog_dist: f64,
    peri_dist: f64,
    focal: bool,
    helio: bool,
    bary: bool,
    earth_center_ecl: [f64; 3],
    earth_bary_ecl: [f64; 3],
    j2000: bool,
    jd_tt: f64,
    ayanamsha_deg: Option<f64>,
) -> Result<[PositionResult; 4]> {
    for v in [node_lon, node_lat, apog_lon, apog_lat, apog_dist, peri_dist] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    if !apog_dist.is_finite() || apog_dist <= 0.0 || !peri_dist.is_finite() || peri_dist <= 0.0 {
        return Err(Error::NonFiniteArgument {
            value: apog_dist.min(peri_dist),
        });
    }
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let peri_lon = (apog_lon + 180.0).rem_euclid(360.0);
    let peri_lat = -apog_lat;
    let a_orb = 0.5 * (peri_dist + apog_dist);
    let e_orb = (apog_dist - peri_dist) / (apog_dist + peri_dist);
    let p_orb = a_orb * (1.0 - e_orb * e_orb);
    let u_node = (node_lon - peri_lon).rem_euclid(360.0).to_radians();
    let nasc_dist = p_orb / (1.0 + e_orb * u_node.cos());
    let ndsc_dist = p_orb / (1.0 - e_orb * u_node.cos());

    let assemble = |lon_d: f64, lat_d: f64, dist: f64| -> Result<PositionResult> {
        let (mut lon, mut lat, mut d) = (lon_d, lat_d, dist);
        if helio || bary {
            // Compatibility contract: FLG_HELCTR adds the heliocentric
            // Earth vector, FLG_BARYCTR the barycentric one; with both
            // set the heliocentric view wins.
            let center = if helio {
                earth_center_ecl
            } else {
                earth_bary_ecl
            };
            let cl = lat.to_radians().cos();
            let gx = d * cl * lon.to_radians().cos() + center[0];
            let gy = d * cl * lon.to_radians().sin() + center[1];
            let gz = d * lat.to_radians().sin() + center[2];
            d = (gx * gx + gy * gy + gz * gz).sqrt();
            lon = gy.atan2(gx).to_degrees().rem_euclid(360.0);
            lat = (gz / d).clamp(-1.0, 1.0).asin().to_degrees();
        }
        if j2000 {
            (lon, lat) = nodaps_dir_to_j2000(lon, lat, jd_tt)?;
        }
        if let Some(aya) = ayanamsha_deg {
            lon = (lon - aya).rem_euclid(360.0);
        }
        Ok([lon, lat, d, 0.0, 0.0, 0.0])
    };

    let xnasc = assemble(node_lon, node_lat, nasc_dist)?;
    let xndsc = assemble((node_lon + 180.0).rem_euclid(360.0), -node_lat, ndsc_dist)?;
    let xperi = assemble(peri_lon, peri_lat, peri_dist)?;
    let xaphe = if focal {
        assemble(apog_lon, apog_lat, apog_dist - peri_dist)?
    } else {
        assemble(apog_lon, apog_lat, apog_dist)?
    };
    Ok([xnasc, xndsc, xperi, xaphe])
}
