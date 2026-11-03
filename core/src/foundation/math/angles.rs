//! Angular normalization and spherical/Cartesian transforms.
//!
//! Ordered layers: radian normalization helpers ([`normalize_2pi`] /
//! [`normalize_pi`]), degree/radian conversion, and spherical/Cartesian
//! transforms both ways over [`crate::math::linalg::Vec3`].
//!
//! Frames, units, epochs: angles are radians unless the name says `deg`;
//! spherical coordinates are `(lon, lat, r)` with `lon = atan2(y, x)` in
//! `[0, 2pi)`, `lat = asin(z / r)` in `[-pi/2, pi/2]`, `r >= 0`. There
//! are no frames or epochs in this module.
//!
//! Conventions (recorded explicitly):
//!
//! - [`normalize_2pi`] maps to `[0, 2pi)`; negative inputs wrap forward and
//!   `2pi` itself maps to `0`. Negative zero maps to `0.0` (positive zero).
//! - [`normalize_pi`] maps to `(-pi, pi]`; `+pi` is kept, `-pi` maps to
//!   `+pi`.
//! - At the poles (`x == 0 && y == 0`, `r > 0`) longitude is defined as `0`
//!   (documented choice, not a geometric fact); the zero vector is an error.
//!
//! Refusals: typed errors, never silent substitution — the zero vector in
//! [`cartesian_to_spherical`], non-finite components or a negative radius
//! in [`spherical_to_cartesian`]. Out-of-range longitudes are accepted and
//! normalized, never refused.
//!
//! Domains and limits: every entry accepts all finite `f64` inputs outside
//! the stated refusals. There is no near-zero snap: small nonzero inputs
//! keep their wrapped value.
//!

use crate::errors::{Error, Result};
use crate::math::linalg::Vec3;
use std::f64::consts::PI;

/// Two pi, the full-circle constant used by the normalization helpers.
pub const TAU: f64 = 2.0 * PI;

/// Degrees per radian.
pub const DEG_PER_RAD: f64 = 180.0 / PI;

/// Normalize an angle in radians to `[0, 2pi)`.
///
/// Uses Euclidean remainder so negative inputs wrap forward; a `2pi`
/// multiple maps to `0`, and negative zero normalizes to positive zero.
/// There is no near-zero snap: small nonzero inputs keep their wrapped value.
pub fn normalize_2pi(mut angle: f64) -> f64 {
    angle = angle.rem_euclid(TAU);
    // rem_euclid can return negative zero for negative-zero input on some
    // paths; canonicalize so `== 0.0` and bit patterns are stable.
    angle + 0.0
}

/// Normalize an angle in radians to `(-pi, pi]`.
///
/// `-pi` maps to `+pi`; this keeps one representative per direction.
pub fn normalize_pi(angle: f64) -> f64 {
    let mut a = normalize_2pi(angle);
    if a > PI {
        a -= TAU;
    }
    // `normalize_2pi` never returns a value in (pi, 2pi) here; the `> PI`
    // branch leaves `(-pi, pi]` with `+pi` preserved and `-pi` mapped up.
    if a == -PI { PI } else { a }
}

/// Radians to degrees.
pub fn rad_to_deg(rad: f64) -> f64 {
    rad * DEG_PER_RAD
}

/// Degrees to radians.
pub fn deg_to_rad(deg: f64) -> f64 {
    deg / DEG_PER_RAD
}

/// Spherical position: longitude in `[0, 2pi)`, latitude in `[-pi/2, pi/2]`,
/// radius `>= 0`. Units of `lon`/`lat` are radians; `r` is the caller's
/// length unit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Spherical {
    /// Longitude `atan2(y, x)`, normalized to `[0, 2pi)`; `0` at the poles.
    pub lon: f64,
    /// Latitude `asin(z / r)`.
    pub lat: f64,
    /// Radius (non-negative).
    pub r: f64,
}

/// Cartesian `(x, y, z)` to spherical. The zero vector is an error; poles
/// yield `lon == 0` by the documented convention above.
pub fn cartesian_to_spherical(v: Vec3) -> Result<Spherical> {
    let (x, y, z) = (v.xyz[0], v.xyz[1], v.xyz[2]);
    let r = (x * x + y * y + z * z).sqrt();
    if r == 0.0 {
        return Err(Error::ZeroVectorToSpherical);
    }
    let lon = if x == 0.0 && y == 0.0 {
        0.0
    } else {
        normalize_2pi(y.atan2(x))
    };
    Ok(Spherical {
        lon,
        lat: (z / r).asin(),
        r,
    })
}

/// Spherical `(lon, lat, r)` to Cartesian. `lon` is normalized first so
/// callers may pass any finite longitude; non-finite inputs are an error.
/// A negative `r` is rejected with [`Error::NonFiniteArgument`]; `lat`
/// outside `[-pi/2, pi/2]` wraps through `sin_cos` without an additional
/// range error.
pub fn spherical_to_cartesian(s: Spherical) -> Result<Vec3> {
    if !s.lon.is_finite() || !s.lat.is_finite() || !s.r.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: s.lon + s.lat + s.r,
        });
    }
    if s.r < 0.0 {
        return Err(Error::NonFiniteArgument { value: s.r });
    }
    let lon = normalize_2pi(s.lon);
    let (sin_lon, cos_lon) = lon.sin_cos();
    let (sin_lat, cos_lat) = s.lat.sin_cos();
    Ok(Vec3::new(
        s.r * cos_lat * cos_lon,
        s.r * cos_lat * sin_lon,
        s.r * sin_lat,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4};

    const TOL: f64 = 1e-12;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= TOL * b.abs().max(1.0)
    }

    #[test]
    fn wrap_boundaries_2pi() {
        assert_eq!(normalize_2pi(0.0), 0.0);
        assert!(close(normalize_2pi(-0.1), TAU - 0.1));
        assert!(close(normalize_2pi(-PI), PI));
        assert_eq!(normalize_2pi(TAU), 0.0);
        assert!(close(normalize_2pi(3.0 * TAU + 1.0), 1.0));
        // Negative zero canonicalizes to positive zero (bit-stable).
        assert_eq!(normalize_2pi(-0.0).to_bits(), 0.0_f64.to_bits());
    }

    #[test]
    fn wrap_boundaries_pi() {
        assert_eq!(normalize_pi(0.0), 0.0);
        assert_eq!(normalize_pi(PI), PI);
        assert_eq!(normalize_pi(-PI), PI);
        assert!(close(normalize_pi(-PI - 0.1), PI - 0.1));
        assert!(close(normalize_pi(PI + 0.1), -PI + 0.1));
        assert!(close(normalize_pi(3.0 * PI), PI));
    }

    #[test]
    fn negative_angles_and_deg_conversions() {
        assert!(close(normalize_2pi(-FRAC_PI_2), 3.0 * FRAC_PI_2));
        assert!(close(rad_to_deg(PI), 180.0));
        assert!(close(deg_to_rad(90.0), FRAC_PI_2));
        assert!(close(rad_to_deg(deg_to_rad(-123.456)), -123.456));
    }

    #[test]
    fn round_trip_off_axis_point() {
        // Non-axis-aligned point: lon = atan2(2, 3), known lat/r checked
        // against direct evaluation, then inverted.
        let v = Vec3::new(3.0, 2.0, -1.0);
        let s = cartesian_to_spherical(v).unwrap();
        assert!(close(s.lon, 2.0_f64.atan2(3.0)));
        let r = 14.0_f64.sqrt();
        assert!(close(s.r, r));
        assert!(close(s.lat, (-1.0 / r).asin()));
        let back = spherical_to_cartesian(s).unwrap();
        assert!(close(back.xyz[0], 3.0));
        assert!(close(back.xyz[1], 2.0));
        assert!(close(back.xyz[2], -1.0));
    }

    #[test]
    fn poles_and_zero_vector() {
        let north = cartesian_to_spherical(Vec3::new(0.0, 0.0, 5.0)).unwrap();
        assert_eq!(north.lon, 0.0);
        assert!(close(north.lat, FRAC_PI_2));
        assert_eq!(north.r, 5.0);
        let south = cartesian_to_spherical(Vec3::new(0.0, -0.0, -2.0)).unwrap();
        assert_eq!(south.lon, 0.0);
        assert!(close(south.lat, -FRAC_PI_2));
        assert!(matches!(
            cartesian_to_spherical(Vec3::zero()),
            Err(Error::ZeroVectorToSpherical)
        ));
    }

    #[test]
    fn spherical_accepts_wrapped_longitude() {
        // lon = -pi/4 is equivalent to 7pi/4: same Cartesian as lon = 7pi/4.
        let a = spherical_to_cartesian(Spherical {
            lon: -FRAC_PI_4,
            lat: 0.0,
            r: 1.0,
        })
        .unwrap();
        let b = spherical_to_cartesian(Spherical {
            lon: TAU - FRAC_PI_4,
            lat: 0.0,
            r: 1.0,
        })
        .unwrap();
        assert!(close(a.xyz[0], b.xyz[0]));
        assert!(close(a.xyz[1], b.xyz[1]));
        assert!(close(a.xyz[2], b.xyz[2]));
        assert!(close(a.xyz[0], FRAC_PI_4.cos()));
        assert!(close(a.xyz[1], -FRAC_PI_4.sin()));
    }
}
