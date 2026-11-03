// SPDX-License-Identifier: AGPL-3.0
//! Maintenance-only osculating elements from a heliocentric J2000 state.
//!
//! The caller owns the TDB epoch and SPK admission. This reducer neither
//! converts clocks nor installs rows in the runtime's element catalogue.

const AU_KM: f64 = 149_597_870.7;
const OBLIQUITY_DEG: f64 = 23.439_291_1;
// Historical Gaussian solar parameter declared by the maintenance artifact,
// not a current nominal solar GM or a mass inferred from the source kernel.
const GM: f64 = 2.959_122_082_855_909_3e-4;

/// Classical elements in the fixed 23.4392911-degree J2000 ecliptic.
///
/// `a` is signed AU, `e` is dimensionless, and `n` is positive degrees/day.
/// Inclination is in [0, 180] degrees; node and perihelion argument are in
/// [0, 360). Elliptic mean anomaly is in [0, 360); hyperbolic mean anomaly
/// is signed and nonperiodic, so reducing it modulo 360 would lose time.
#[derive(Clone, Copy, Debug)]
pub(super) struct Elements {
    pub a: f64,
    pub e: f64,
    pub i: f64,
    pub omega: f64,
    pub node: f64,
    pub m0: f64,
    pub n: f64,
}

/// Reduce finite Sun-relative J2000-equatorial km and km/s to elements.
///
/// Both vectors receive the same passive x rotation. Geometry follows
/// h = r × v, e = (v × h)/GM − r/|r| and a = −GM/(2 energy).
/// Oriented angles use atan2 in the orbital plane, retaining retrograde
/// orientation. For an exactly equatorial orbit the reference node is +x; at
/// exactly zero eccentricity perihelion is that node and mean anomaly is
/// the argument of latitude. Small nonzero inclinations/eccentricities are
/// not rounded away by an arbitrary physical threshold.
///
/// Zero radius/angular momentum, parabolic energy, inconsistent unresolved
/// conic classification, and nonrepresentable intermediates are errors.
/// This is a binary64 maintenance calculation, not a universal conditioning
/// guarantee for states arbitrarily close to radial/parabolic degeneracies.
pub(super) fn from_equatorial_state(
    position_km: [f64; 3],
    velocity_km_s: [f64; 3],
) -> Result<Elements, String> {
    if !position_km
        .into_iter()
        .chain(velocity_km_s)
        .all(f64::is_finite)
    {
        return Err("nonfinite Cartesian state".to_owned());
    }
    let (sin_eps, cos_eps) = OBLIQUITY_DEG.to_radians().sin_cos();
    let rotate = |v: [f64; 3]| {
        [
            v[0],
            v[1] * cos_eps + v[2] * sin_eps,
            -v[1] * sin_eps + v[2] * cos_eps,
        ]
    };
    from_ecliptic_state(
        rotate(position_km.map(|value| value / AU_KM)),
        rotate(velocity_km_s.map(|value| value / AU_KM * 86_400.0)),
    )
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(v: [f64; 3]) -> f64 {
    v[0].hypot(v[1]).hypot(v[2])
}

fn angle(a: [f64; 3], b: [f64; 3], normal: [f64; 3]) -> f64 {
    dot(cross(a, b), normal).atan2(dot(a, b))
}

fn degrees(angle: f64) -> f64 {
    let value = angle.to_degrees().rem_euclid(360.0);
    if value == 360.0 || value == 0.0 {
        0.0
    } else {
        value
    }
}

fn from_ecliptic_state(r: [f64; 3], v: [f64; 3]) -> Result<Elements, String> {
    let radius = norm(r);
    let h = cross(r, v);
    let angular_momentum = norm(h);
    let energy = dot(v, v) / 2.0 - GM / radius;
    if !radius.is_finite()
        || radius == 0.0
        || !angular_momentum.is_finite()
        || angular_momentum == 0.0
        || !energy.is_finite()
    {
        return Err("degenerate or nonrepresentable Cartesian state".to_owned());
    }
    if energy == 0.0 {
        return Err("parabolic states have no finite semimajor axis".to_owned());
    }
    let a = -GM / (2.0 * energy);
    let r_unit = r.map(|value| value / radius);
    let h_unit = h.map(|value| value / angular_momentum);
    let vh = cross(v, h);
    let eccentricity = std::array::from_fn(|k| vh[k] / GM - r_unit[k]);
    let e = norm(eccentricity);
    if !a.is_finite()
        || a == 0.0
        || !e.is_finite()
        || (energy < 0.0 && e >= 1.0)
        || (energy > 0.0 && e <= 1.0)
    {
        return Err("nonrepresentable or numerically unresolved conic".to_owned());
    }
    let node_norm = h[0].hypot(h[1]);
    let node_unit = if node_norm == 0.0 {
        [1.0, 0.0, 0.0]
    } else {
        [-h[1] / node_norm, h[0] / node_norm, 0.0]
    };
    let (omega, mean_anomaly) = if e == 0.0 {
        (0.0, angle(node_unit, r_unit, h_unit))
    } else {
        let perihelion = eccentricity.map(|value| value / e);
        let omega = angle(node_unit, perihelion, h_unit);
        let mean = if e < 1.0 {
            let nu = angle(perihelion, r_unit, h_unit);
            let (sin_nu, cos_nu) = nu.sin_cos();
            let anomaly = (((1.0 - e) * (1.0 + e)).sqrt() * sin_nu).atan2(e + cos_nu);
            (-e).mul_add(anomaly.sin(), anomaly)
        } else {
            // r·v = e sqrt(GM |a|) sinh(H); no asymptote-angle division.
            let sinh_h = (dot(r, v) / (GM * -a).sqrt()) / e;
            e.mul_add(sinh_h, -sinh_h.asinh())
        };
        (omega, mean)
    };
    let result = Elements {
        a,
        e,
        i: node_norm.atan2(h[2]).to_degrees(),
        omega: degrees(omega),
        node: if node_norm == 0.0 {
            0.0
        } else {
            degrees(h[0].atan2(-h[1]))
        },
        m0: if e < 1.0 {
            degrees(mean_anomaly)
        } else {
            mean_anomaly.to_degrees()
        },
        n: (GM.sqrt() / a.abs().sqrt() / a.abs()).to_degrees(),
    };
    if ![
        result.a,
        result.e,
        result.i,
        result.omega,
        result.node,
        result.m0,
        result.n,
    ]
    .into_iter()
    .all(f64::is_finite)
        || result.n <= 0.0
    {
        return Err("nonrepresentable orbital elements".to_owned());
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Independent forward Kepler equations in the perifocal plane, followed
    // by the three Euler rotations; no Cartesian-to-elements helper is used.
    fn state(a: f64, e: f64, i: f64, node: f64, omega: f64, anomaly: f64) -> ([f64; 3], [f64; 3]) {
        let (x, y, vx, vy) = if e < 1.0 {
            let (sin_e, cos_e) = anomaly.sin_cos();
            let beta = ((1.0 - e) * (1.0 + e)).sqrt();
            let rate = (GM / a).sqrt() / (1.0 - e * cos_e);
            (
                a * (cos_e - e),
                a * beta * sin_e,
                -rate * sin_e,
                rate * beta * cos_e,
            )
        } else {
            let (sin_h, cos_h) = (anomaly.sinh(), anomaly.cosh());
            let beta = ((e - 1.0) * (e + 1.0)).sqrt();
            let rate = (GM / -a).sqrt() / (e * cos_h - 1.0);
            (
                -a * (e - cos_h),
                -a * beta * sin_h,
                -rate * sin_h,
                rate * beta * cos_h,
            )
        };
        let (sn, cn) = node.sin_cos();
        let (si, ci) = i.sin_cos();
        let (sw, cw) = omega.sin_cos();
        let rotate = |x: f64, y: f64| {
            let p = cw * x - sw * y;
            let q = sw * x + cw * y;
            [cn * p - sn * ci * q, sn * p + cn * ci * q, si * q]
        };
        (rotate(x, y), rotate(vx, vy))
    }

    fn equatorial(r: [f64; 3], v: [f64; 3]) -> ([f64; 3], [f64; 3]) {
        let (s, c) = OBLIQUITY_DEG.to_radians().sin_cos();
        let rotate = |p: [f64; 3]| [p[0], c * p[1] - s * p[2], s * p[1] + c * p[2]];
        (
            rotate(r).map(|x| x * AU_KM),
            rotate(v).map(|x| x * AU_KM / 86_400.0),
        )
    }

    #[test]
    fn known_perihelion_retains_scale_mass_and_ecliptic_convention() {
        let r = [1.5, 0.0, 0.0];
        let v = [0.0, (GM * (2.0 / 1.5 - 1.0 / 2.0)).sqrt(), 0.0];
        let (r_km, v_km_s) = equatorial(r, v);
        let elements = from_equatorial_state(r_km, v_km_s).unwrap();
        assert!((elements.a - 2.0).abs() < 2e-14);
        assert!((elements.e - 0.25).abs() < 2e-14);
        assert!(elements.i.abs() < 2e-13);
        assert!((elements.n - (GM / 8.0).sqrt().to_degrees()).abs() < 2e-15);
        assert!(elements.m0.min(360.0 - elements.m0) < 2e-12);
    }

    #[test]
    fn reconstructed_states_cover_circular_inclined_retrograde_and_hyperbolic_orbits() {
        for (a, e, i, node, omega, anomaly) in [
            (2.5, 0.4, 0.7, 1.2, 0.8, -1.3),
            (4.0, 0.8, 2.7, 2.4, 1.1, 2.1),
            (1.0, 0.0, 0.0, 0.0, 0.0, 0.7),
            (1.0, 0.0, 1.2, 0.5, 0.0, -0.4),
            (2.0, 0.2, 0.0, 0.0, 0.9, 0.8),
            (2.0, 0.2, std::f64::consts::PI, 0.0, 0.9, 0.8),
            (-3.0, 1.6, 0.9, 0.2, 2.0, -1.5),
            (-3.0, 1.6, 2.4, 0.2, 2.0, 3.0),
        ] {
            let (r, v) = state(a, e, i, node, omega, anomaly);
            let (r_km, v_km_s) = equatorial(r, v);
            let el = from_equatorial_state(r_km, v_km_s).unwrap();
            assert!((el.a - a).abs() / a.abs() < 2e-13, "{el:?}");
            assert!((el.e - e).abs() < 2e-13, "{el:?}");
            let mean = el.m0.to_radians();
            let mut recovered = if el.e < 1.0 {
                mean
            } else {
                (mean / el.e).asinh()
            };
            for _ in 0..32 {
                recovered -= if el.e < 1.0 {
                    (recovered - el.e * recovered.sin() - mean) / (1.0 - el.e * recovered.cos())
                } else {
                    (el.e * recovered.sinh() - recovered - mean) / (el.e * recovered.cosh() - 1.0)
                };
            }
            let (rr, vv) = state(
                el.a,
                el.e,
                el.i.to_radians(),
                el.node.to_radians(),
                el.omega.to_radians(),
                recovered,
            );
            assert!(
                norm(std::array::from_fn(|k| rr[k] - r[k])) / norm(r) < 2e-13,
                "{el:?}"
            );
            assert!(
                norm(std::array::from_fn(|k| vv[k] - v[k])) / norm(v) < 2e-13,
                "{el:?}"
            );
            if e > 1.0 {
                assert_eq!(el.m0.is_sign_negative(), anomaly.is_sign_negative());
            }
        }
    }

    #[test]
    fn invalid_radial_and_exact_parabolic_states_are_refused() {
        assert!(from_equatorial_state([0.0; 3], [1.0; 3]).is_err());
        assert!(from_equatorial_state([AU_KM, 0.0, 0.0], [1.0, 0.0, 0.0]).is_err());
        assert!(from_ecliptic_state([2.0 * GM, 0.0, 0.0], [0.0, 1.0, 0.0]).is_err());
        for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(from_equatorial_state([invalid, 0.0, 0.0], [0.0, 1.0, 0.0]).is_err());
            assert!(from_equatorial_state([AU_KM, 0.0, 0.0], [0.0, invalid, 0.0]).is_err());
        }
        assert!(from_equatorial_state([f64::MAX; 3], [f64::MAX; 3]).is_err());
    }
}
