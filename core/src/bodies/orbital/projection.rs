//! Node/apse frame rotations and observer projection.
//!
//! Projection uses borrowed deflector callbacks only under the established
//! physical gates. It does not acquire targets or choose serving clocks.

use super::{
    J2000_JD, NODAPS_SPEED_STEP_DAYS, NOMINAL_SOLAR_RADIUS_AU, NodapsDeflectorAt, NodapsFrame,
    NodapsProjection,
};
use crate::constants;
use crate::errors::{Error, Result};
use crate::positions::{
    C_LIGHT_AU_DAY, PositionResult, apply_aberration, mat3_vec3, vondrak_mean_obliquity_rad,
    vondrak_precession_matrix,
};

impl NodapsFrame {
    /// Build the working ecliptic frame at TT Julian day `jd_tt`.
    /// `nutation` supplies longitude/obliquity in radians; `None` selects mean
    /// coordinates. Nonfinite epochs or frame components return a domain error.
    pub fn for_date(jd_tt: f64, nutation: Option<(f64, f64)>) -> Result<Self> {
        let (dpsi, deps) = nutation.unwrap_or((0.0, 0.0));
        let frames = crate::positions::frame_matrices_for_date(jd_tt, dpsi, deps)?;
        let flags = if nutation.is_some() {
            0
        } else {
            constants::FLG_NONUT
        };
        let mut columns = [[0.0; 3]; 3];
        for (column, output) in columns.iter_mut().enumerate() {
            let mut basis = [0.0; 3];
            basis[column] = 1.0;
            *output = crate::positions::frame_transform(basis, jd_tt, flags, true, frames)?;
        }
        Ok(Self {
            matrix: std::array::from_fn(|row| std::array::from_fn(|column| columns[column][row])),
        })
    }

    /// Rotate ICRS components into ecliptic-of-date components, preserving units.
    /// Input components share arbitrary units; finite input yields a finite
    /// rotation result unless arithmetic overflows. No data are acquired.
    pub fn icrs_to_ecliptic(self, vector: [f64; 3]) -> [f64; 3] {
        self.matrix
            .map(|row| crate::positions::compensated_dot3(row, vector))
    }

    /// Rotate ecliptic-of-date components back to ICRS, preserving units.
    /// The transpose is the inverse of the orthonormal physical frame, up to
    /// its binary64 construction error. Nonfinite input propagates naturally.
    pub fn ecliptic_to_icrs(self, vector: [f64; 3]) -> [f64; 3] {
        std::array::from_fn(|column| {
            crate::positions::compensated_dot3(
                std::array::from_fn(|row| self.matrix[row][column]),
                vector,
            )
        })
    }
}

/// Rotate an ecliptic-of-date vector into the node/apse J2000 frame
/// : to equatorial of
/// date with the J2000 mean obliquity, through the transposed biased
/// date precession, through the biased J2000 precession, and back to
/// the J2000 ecliptic. Linear, so position and velocity vectors
/// transform identically. Errors: returns `Error::NonFiniteArgument` on
/// any non-finite component or `jd_tt` (a checked-construction addition;
/// the model chain itself is total on finite inputs). Propagates through
/// `nodaps_dir_to_j2000`.
pub fn nodaps_vector_to_j2000(v: [f64; 3], jd_tt: f64) -> Result<[f64; 3]> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    for c in v {
        if !c.is_finite() {
            return Err(Error::NonFiniteArgument { value: c });
        }
    }
    let eps_j = vondrak_mean_obliquity_rad(J2000_JD)?;
    let p_date = vondrak_precession_matrix(jd_tt, true)?;
    let p_j2k = vondrak_precession_matrix(J2000_JD, true)?;
    let (sin_j, cos_j) = eps_j.sin_cos();
    // Ecliptic-of-date -> equatorial-of-date (J2000 obliquity).
    let mut w = [
        v[0],
        v[1] * cos_j - v[2] * sin_j,
        v[1] * sin_j + v[2] * cos_j,
    ];
    // Equator of date -> ICRS (transpose) -> equator J2000.
    w = [
        p_date[0][0] * w[0] + p_date[1][0] * w[1] + p_date[2][0] * w[2],
        p_date[0][1] * w[0] + p_date[1][1] * w[1] + p_date[2][1] * w[2],
        p_date[0][2] * w[0] + p_date[1][2] * w[1] + p_date[2][2] * w[2],
    ];
    w = [
        p_j2k[0][0] * w[0] + p_j2k[0][1] * w[1] + p_j2k[0][2] * w[2],
        p_j2k[1][0] * w[0] + p_j2k[1][1] * w[1] + p_j2k[1][2] * w[2],
        p_j2k[2][0] * w[0] + p_j2k[2][1] * w[1] + p_j2k[2][2] * w[2],
    ];
    // Equatorial -> ecliptic (J2000 obliquity).
    Ok([
        w[0],
        w[1] * cos_j + w[2] * sin_j,
        -w[1] * sin_j + w[2] * cos_j,
    ])
}

/// Precess a node/apse of-date ecliptic direction to the J2000 frame
/// : the vector chain above, then
/// longitude/latitude extraction. Nutation present in the of-date
/// coordinate is deliberately retained. Returns `(lon_deg, lat_deg)`.
pub fn nodaps_dir_to_j2000(lon_deg: f64, lat_deg: f64, jd_tt: f64) -> Result<(f64, f64)> {
    let cl = lat_deg.to_radians().cos();
    let v = nodaps_vector_to_j2000(
        [
            cl * lon_deg.to_radians().cos(),
            cl * lon_deg.to_radians().sin(),
            lat_deg.to_radians().sin(),
        ],
        jd_tt,
    )?;
    let lon_j = v[1].atan2(v[0]).to_degrees().rem_euclid(360.0);
    let r = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    let lat_j = (v[2] / r).clamp(-1.0, 1.0).asin().to_degrees();
    Ok((lon_j, lat_j))
}

/// Whether the finite observer-to-target ray intersects the nominal
/// solar disk: inside the
/// disk the point-mass PPN reduction is physically invalid, so the
/// deflection is skipped for that ray.
pub fn nodaps_ray_crosses_solar_disk(
    observer_to_target: [f64; 3],
    observer_barycentric: [f64; 3],
    sun_barycentric: [f64; 3],
) -> bool {
    let distance = (observer_to_target[0] * observer_to_target[0]
        + observer_to_target[1] * observer_to_target[1]
        + observer_to_target[2] * observer_to_target[2])
        .sqrt();
    if distance == 0.0 {
        return false;
    }
    let direction = [
        observer_to_target[0] / distance,
        observer_to_target[1] / distance,
        observer_to_target[2] / distance,
    ];
    let observer_to_sun = [
        sun_barycentric[0] - observer_barycentric[0],
        sun_barycentric[1] - observer_barycentric[1],
        sun_barycentric[2] - observer_barycentric[2],
    ];
    let along_ray = observer_to_sun[0] * direction[0]
        + observer_to_sun[1] * direction[1]
        + observer_to_sun[2] * direction[2];
    if along_ray <= 0.0 || along_ray >= distance {
        return false;
    }
    let impact = ((observer_to_sun[0] - along_ray * direction[0]).powi(2)
        + (observer_to_sun[1] - along_ray * direction[1]).powi(2)
        + (observer_to_sun[2] - along_ray * direction[2]).powi(2))
    .sqrt();
    impact < NOMINAL_SOLAR_RADIUS_AU
}

impl<'d> NodapsProjection<'d> {
    /// Checked constructor. `observer_ecl` is the vector subtracted from
    /// the orbit-frame point, or `None` for a point already relative to
    /// its requested origin. A heliocentric fit viewed from the SSB uses
    /// the negative barycentric Sun as this offset. `mirror_sun` enables
    /// the solar-orbit negation only for a non-centered solar view; the
    /// caller must clear it for heliocentric and barycentric requests.
    /// `mean_align_dpsi` carries the
    /// nutation-in-longitude for the mean-element frame alignment
    /// (`None` under `FLG_NONUT` or for osculating shapes, which need no
    /// alignment). Errors: returns `Error::NonFiniteArgument` on any
    /// non-finite observer component, earth-velocity/barycentric component,
    /// sun-barycentric component, `jd_tt`, or `mean_align_dpsi` value.
    /// Positions and projected points are AU; `earth_vel_ecl` is AU/day in the
    /// working ecliptic frame. Barycentric Earth/Sun positions are ICRS AU.
    /// `jd_tt` is a TT Julian day and `mean_align_dpsi` is radians.
    /// `deflect_at(body_id, jd_tt)` returns the deflector's barycentric ICRS
    /// position in AU at each requested TT epoch. The callback is borrowed for
    /// this projection's lifetime and runs only if a point consumes deflection;
    /// consumed callback errors propagate without replacement states.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        observer_ecl: Option<[f64; 3]>,
        mirror_sun: bool,
        apply_aberration: bool,
        apply_deflection: bool,
        earth_vel_ecl: [f64; 3],
        earth_bary_icrs: [f64; 3],
        sun_bary_icrs: [f64; 3],
        jd_tt: f64,
        frame: NodapsFrame,
        mean_align_dpsi: Option<f64>,
        deflect_at: &'d dyn Fn(i32, f64) -> Result<[f64; 3]>,
    ) -> Result<Self> {
        if let Some(obs) = observer_ecl {
            for c in obs {
                if !c.is_finite() {
                    return Err(Error::NonFiniteArgument { value: c });
                }
            }
        }
        for c in earth_vel_ecl
            .into_iter()
            .chain(earth_bary_icrs)
            .chain(sun_bary_icrs)
        {
            if !c.is_finite() {
                return Err(Error::NonFiniteArgument { value: c });
            }
        }
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        if let Some(dpsi) = mean_align_dpsi {
            if !dpsi.is_finite() {
                return Err(Error::NonFiniteArgument { value: dpsi });
            }
        }
        Ok(Self {
            observer_ecl,
            mirror_sun,
            apply_aberration,
            apply_deflection,
            earth_vel_ecl,
            earth_bary_icrs,
            sun_bary_icrs,
            jd_tt,
            frame,
            mean_align_dpsi,
            deflect_at,
            optional_deflect_at: None,
        })
    }

    /// Attach a source adapter whose raw missing/range deflectors are optional.
    /// The callback receives scalar TT observation/closest-approach epochs and
    /// returns an SSB ICRS AU vector, or None for a source-skipped deflector.
    /// Other typed errors propagate. The numerical PPN kernel and default
    /// strict lookup remain unchanged; the borrowed feed cannot outlive this
    /// projection and performs no apparent target acquisition.
    pub(super) fn with_optional_deflector_lookup(
        mut self,
        lookup: &'d NodapsDeflectorAt<'d>,
    ) -> Self {
        self.optional_deflect_at = Some(lookup);
        self
    }

    /// Project one orbit-frame point.
    /// Returns `(lon_deg, lat_deg, dist_au)`. Branches: the mean-alignment
    /// rotation applies first when set; the mirrored Sun negates only for
    /// non-centered requests (centered requests keep the orbit point);
    /// deflection is skipped when the ray crosses the solar disk; a
    /// zero-distance point yields latitude 0.0.
    pub fn project(&self, point: [f64; 3]) -> Result<(f64, f64, f64)> {
        let mut p = point;
        if let Some(dpsi) = self.mean_align_dpsi {
            let (sin_d, cos_d) = dpsi.sin_cos();
            p = [
                p[0] * cos_d - p[1] * sin_d,
                p[0] * sin_d + p[1] * cos_d,
                p[2],
            ];
        }
        let (mut gx, mut gy, mut gz) = if self.mirror_sun && self.observer_ecl.is_some() {
            // Apparent solar orbit: the heliocentric Earth-orbit point
            // mirrored through the origin. The caller has suppressed
            // this mirror for every centered request, independently of
            // whether that view needs an observer offset.
            (-p[0], -p[1], -p[2])
        } else if let Some(obs) = self.observer_ecl {
            (p[0] - obs[0], p[1] - obs[1], p[2] - obs[2])
        } else {
            (p[0], p[1], p[2])
        };
        if self.apply_deflection {
            let geo_icrs = self.frame.ecliptic_to_icrs([gx, gy, gz]);
            if !nodaps_ray_crosses_solar_disk(geo_icrs, self.earth_bary_icrs, self.sun_bary_icrs) {
                let distance = (gx * gx + gy * gy + gz * gz).sqrt();
                let bent = crate::positions::apply_gravitational_deflection_with_lookup(
                    geo_icrs,
                    self.earth_bary_icrs,
                    self.jd_tt,
                    distance / C_LIGHT_AU_DAY,
                    |id, jd| match self.optional_deflect_at {
                        Some(lookup) => lookup(id, jd),
                        None => (self.deflect_at)(id, jd).map(Some),
                    },
                )?;
                let back = self.frame.icrs_to_ecliptic(bent);
                gx = back[0];
                gy = back[1];
                gz = back[2];
            }
        }
        if self.apply_aberration {
            // `apply_aberration_to_position` (first-order Bradley
            // form): the native with zero light time is its
            // structural equivalent (floors differ only below 1 e-30 AU,
            // unphysical for node/apse points at AU scales).
            let bent = apply_aberration([gx, gy, gz], self.earth_vel_ecl, 0.0)?;
            gx = bent[0];
            gy = bent[1];
            gz = bent[2];
        }
        let r_geo = (gx * gx + gy * gy + gz * gz).sqrt();
        let lon = gy.atan2(gx).to_degrees().rem_euclid(360.0);
        let lat = if r_geo > 0.0 {
            (gz / r_geo).clamp(-1.0, 1.0).asin().to_degrees()
        } else {
            0.0
        };
        Ok((lon, lat, r_geo))
    }

    /// Project acquired conic points in ascending-node, descending-node,
    /// perihelion and aphelion order. The solar/terrestrial orbit defines the
    /// ecliptic: its undefined nodes return zeros without invoking projection
    /// or a deflector feed. Callback composition and later output-frame or
    /// ayanamsha reductions remain the serving layer's responsibility.
    pub(super) fn project_orbit_points(
        &self,
        points: [[f64; 3]; 4],
        body_id: i32,
        mirror_sun: bool,
    ) -> Result<[(f64, f64, f64); 4]> {
        let [asc, dsc, peri, aphe] = points;
        let (asc, dsc) = if body_id == constants::EARTH || mirror_sun {
            ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
        } else {
            (self.project(asc)?, self.project(dsc)?)
        };
        Ok([asc, dsc, self.project(peri)?, self.project(aphe)?])
    }
}

/// FLG_J2000 longitude rule for finished of-date geocentric points
/// : the reduction precesses only
/// the ecliptic longitude of each heliocentric node/apse by the
/// accumulated general precession, keeping the of-date latitude, then
/// re-projects geocentrically against the J2000 observer. Points with
/// zero distance (Sun/Earth node sentinels) pass through untouched.
/// `observer_ecl` is the same observer subtracted above (`None` for
/// centered requests, where the correction vector vanishes, as for the
/// mirrored Sun). The frame-to-frame rotation is mean-ecliptic only
/// (accumulated general precession, nutation excluded); the of-date
/// coordinates themselves retain their nutation. Re-projection against the
/// J2000 observer leaves a distance-dependent parallactic latitude (near
/// zero at J2000, growing with `|t - J2000|`) that a rigid rotation of the
/// geocentric direction cannot reproduce. Errors: propagates
/// frame-construction failures on non-finite `jd_tt`.
pub fn apply_nodaps_j2000_rule(
    points: [(f64, f64, f64); 4],
    observer_ecl: Option<[f64; 3]>,
    mirror_sun: bool,
    jd_tt: f64,
) -> Result<[(f64, f64, f64); 4]> {
    let pn_mean = vondrak_precession_matrix(jd_tt, true)?;
    let eps_mean = vondrak_mean_obliquity_rad(jd_tt)?;
    let (sin_m, cos_m) = eps_mean.sin_cos();
    // Mean-ecliptic-of-date -> ICRS through the biased date precession.
    let mean_ecl_to_icrs = |v: [f64; 3]| -> [f64; 3] {
        let xq = v[0];
        let yq = v[1] * cos_m - v[2] * sin_m;
        let zq = v[1] * sin_m + v[2] * cos_m;
        [
            pn_mean[0][0] * xq + pn_mean[1][0] * yq + pn_mean[2][0] * zq,
            pn_mean[0][1] * xq + pn_mean[1][1] * yq + pn_mean[2][1] * zq,
            pn_mean[0][2] * xq + pn_mean[1][2] * yq + pn_mean[2][2] * zq,
        ]
    };
    // ICRS -> mean ecliptic of J2000.
    let to_ecl_j2000 = |v: [f64; 3]| -> Result<[f64; 3]> {
        let pj = vondrak_precession_matrix(J2000_JD, true)?;
        let eps_j = vondrak_mean_obliquity_rad(J2000_JD)?;
        let (sin_j, cos_j) = eps_j.sin_cos();
        let eq = mat3_vec3(pj, v)?;
        Ok([
            eq[0],
            eq[1] * cos_j + eq[2] * sin_j,
            -eq[1] * sin_j + eq[2] * cos_j,
        ])
    };
    let eqx = to_ecl_j2000(mean_ecl_to_icrs([1.0, 0.0, 0.0]))?;
    let delta = eqx[1].atan2(eqx[0]);
    let (sin_d, cos_d) = delta.sin_cos();
    let cvec = if observer_ecl.is_none() || mirror_sun {
        [0.0; 3]
    } else {
        let obs = observer_ecl.unwrap_or([0.0; 3]);
        let obs_j = to_ecl_j2000(mean_ecl_to_icrs(obs))?;
        [
            cos_d * obs[0] - sin_d * obs[1] - obs_j[0],
            sin_d * obs[0] + cos_d * obs[1] - obs_j[1],
            obs[2] - obs_j[2],
        ]
    };
    let mut out = [(0.0, 0.0, 0.0); 4];
    for (i, g) in points.iter().enumerate() {
        if g.2 == 0.0 {
            out[i] = *g;
            continue;
        }
        let cl = g.1.to_radians().cos();
        let vx = g.2 * cl * g.0.to_radians().cos();
        let vy = g.2 * cl * g.0.to_radians().sin();
        let vz = g.2 * g.1.to_radians().sin();
        let rx = cos_d * vx - sin_d * vy + cvec[0];
        let ry = sin_d * vx + cos_d * vy + cvec[1];
        let rz = vz + cvec[2];
        let rr = (rx * rx + ry * ry + rz * rz).sqrt();
        out[i] = (
            ry.atan2(rx).to_degrees().rem_euclid(360.0),
            (rz / rr).clamp(-1.0, 1.0).asin().to_degrees(),
            rr,
        );
    }
    Ok(out)
}

/// Fourth-order centered derivative of one spherical node or apsis.
///
/// `cur`, `pre`/`nex` and `pre2`/`nex2` are the position at TT, at TT ± h
/// and at TT ± 2 h, where h is [`NODAPS_SPEED_STEP_DAYS`]. Longitudes and
/// latitudes are degrees, distances AU; returned rates are degrees/day and
/// AU/day. Longitude increments take the short arc, preserving a continuous
/// derivative across the zero-degree seam. Other coordinates are ordinary
/// binary64 differences. Nonfinite coordinates propagate.
pub fn nodaps_point_speeds(
    cur: (f64, f64, f64),
    pre: (f64, f64, f64),
    nex: (f64, f64, f64),
    pre2: (f64, f64, f64),
    nex2: (f64, f64, f64),
) -> PositionResult {
    let longitude_delta = |next: f64, previous: f64| {
        let delta = next - previous;
        if delta > 180.0 {
            delta - 360.0
        } else if delta < -180.0 {
            delta + 360.0
        } else {
            delta
        }
    };
    let derivative = |near: f64, far: f64| (8.0 * near - far) / (12.0 * NODAPS_SPEED_STEP_DAYS);
    [
        cur.0,
        cur.1,
        cur.2,
        derivative(
            longitude_delta(nex.0, pre.0),
            longitude_delta(nex2.0, pre2.0),
        ),
        derivative(nex.1 - pre.1, nex2.1 - pre2.1),
        derivative(nex.2 - pre.2, nex2.2 - pre2.2),
    ]
}
