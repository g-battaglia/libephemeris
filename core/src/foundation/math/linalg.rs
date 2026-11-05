//! Fixed-size three-vector and 3x3 matrix operations.
//!
//! Ordered layers: [`Vec3`] (dot/cross products, norms, normalization,
//! component-wise arithmetic) and [`Mat3`] (active right-handed rotations,
//! transpose, product, application, adjugate inverse, orthogonality
//! checks) in the written association order.
//!
//! Frames, units, epochs: vector components and matrix entries carry the
//! caller's meaning and units; matrices are row-major `m[row][col]` with
//! `apply` computing `y_i = sum_j m[i][j] * v[j]`. There are no frames or
//! epochs in this module.
//!
//! These cover what the astrometry pipeline needs: dot/cross products,
//! norms, rotations, transposes, inverses and orthogonality checks. No
//! generic tensor framework: three-vectors and 3x3 matrices are sufficient.
//!
//! Conventions (recorded explicitly; verified by the tests below):
//!
//! - Right-handed rotations about the indexed axis, active convention:
//!   `rot_z(90deg)` maps `(1, 0, 0)` to `(0, 1, 0)`.
//! - Matrices are row-major `m[row][col]`; `apply` computes `y_i =
//!  sum_j m[i][j] * v[j]` in the written association order.
//! - [`Mat3::inverse`] uses the adjugate/determinant form and rejects
//!   singular inputs with [`Error::SingularMatrix`] instead of producing
//!   infinities.
//!
//! Refusals: typed errors, never silent substitution — exact-zero vector
//! normalization and exactly singular (`det == 0.0`) matrix inversion.
//! An overflowing norm scales by `1.0 / inf` and a NaN norm propagates NaN
//! components; callers with extreme magnitudes must check finiteness before
//! normalizing.
//!
//! Domains and limits: every entry accepts all finite `f64` inputs outside
//! the stated refusals; products and determinants accumulate in the written
//! association order.
//!

use crate::errors::{Error, Result};
use std::ops::{Add, Mul, Sub};

/// Euclidean two-component norm with compensated squares and root correction.
///
/// Components share caller-defined units; the nonnegative result has those
/// units. Order absolute magnitudes and divide by a binary power to bound the
/// squares. FMA recovers each product residue, FastTwoSum recovers the addition
/// residue, and `r += (sum - r² + low)/(2r)` corrects the square root before
/// rescaling. This retains low words lost by a platform `hypot` primitive.
/// Zero returns positive zero; infinity takes precedence over NaN. A subnormal
/// maximum scales by the smallest normal power. A true result above the finite
/// binary64 range overflows to infinity. This floating-point primitive supplies
/// no arbitrary-precision proof or universal correctly-rounded guarantee.
pub(crate) fn compensated_hypot(x: f64, y: f64) -> f64 {
    if x.is_infinite() || y.is_infinite() {
        return f64::INFINITY;
    }
    if x.is_nan() || y.is_nan() {
        return f64::NAN;
    }
    let (large, small) = (x.abs().max(y.abs()), x.abs().min(y.abs()));
    if large == 0.0 {
        return 0.0;
    }
    let exponent = large.to_bits() & 0x7ff0_0000_0000_0000;
    let scale = f64::from_bits(exponent.max(0x0010_0000_0000_0000));
    let (a, b) = (large / scale, small / scale);
    let (aa, bb) = (a * a, b * b);
    let sum = aa + bb;
    let low = (aa - sum) + bb + a.mul_add(a, -aa) + b.mul_add(b, -bb);
    let root = sum.sqrt();
    let correction = ((-root).mul_add(root, sum) + low) / (2.0 * root);
    (root + correction) * scale
}

/// Three-component Cartesian vector of `f64`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vec3 {
    /// Cartesian components `(x, y, z)`; meaning/units are the caller's.
    pub xyz: [f64; 3],
}

impl Vec3 {
    /// Construct from components (no validation; meaning is caller-defined).
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { xyz: [x, y, z] }
    }

    /// Zero vector.
    pub fn zero() -> Self {
        Self { xyz: [0.0; 3] }
    }

    /// Scalar multiplication.
    pub fn scale(self, s: f64) -> Self {
        Self {
            xyz: [self.xyz[0] * s, self.xyz[1] * s, self.xyz[2] * s],
        }
    }

    /// Dot product, accumulated in component order.
    pub fn dot(self, other: Self) -> f64 {
        self.xyz[0] * other.xyz[0] + self.xyz[1] * other.xyz[1] + self.xyz[2] * other.xyz[2]
    }

    /// Cross product `self x other` (right-handed).
    pub fn cross(self, other: Self) -> Self {
        let (a, b) = (self.xyz, other.xyz);
        Self {
            xyz: [
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0],
            ],
        }
    }

    /// Euclidean norm.
    pub fn norm(self) -> f64 {
        self.dot(self).sqrt()
    }

    /// Unit vector; rejects the zero vector instead of dividing by zero.
    ///
    /// Only the exact zero vector is rejected: an overflowing norm
    /// (infinite for extreme magnitudes) scales by `1.0 / inf = 0.0` and
    /// returns the zero vector without an error, and a NaN norm propagates
    /// NaN components. Callers with extreme magnitudes must check finiteness
    /// before normalizing.
    pub fn normalized(self) -> Result<Self> {
        let n = self.norm();
        if n == 0.0 {
            return Err(Error::ZeroVectorNormalize);
        }
        Ok(self.scale(1.0 / n))
    }
}

impl Add for Vec3 {
    type Output = Self;

    /// Component-wise sum, in the written association order.
    fn add(self, other: Self) -> Self {
        Self {
            xyz: [
                self.xyz[0] + other.xyz[0],
                self.xyz[1] + other.xyz[1],
                self.xyz[2] + other.xyz[2],
            ],
        }
    }
}

impl Sub for Vec3 {
    type Output = Self;

    /// Component-wise difference.
    fn sub(self, other: Self) -> Self {
        Self {
            xyz: [
                self.xyz[0] - other.xyz[0],
                self.xyz[1] - other.xyz[1],
                self.xyz[2] - other.xyz[2],
            ],
        }
    }
}

/// Row-major 3x3 matrix of `f64`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Mat3 {
    /// Entries `m[row][col]`.
    pub m: [[f64; 3]; 3],
}

impl Mul for Mat3 {
    type Output = Self;

    /// Matrix product, inner accumulation in `k` order.
    fn mul(self, other: Self) -> Self {
        let mut m = [[0.0; 3]; 3];
        for (out_row, self_row) in m.iter_mut().zip(self.m.iter()) {
            for (j, v) in out_row.iter_mut().enumerate() {
                *v = self_row[0] * other.m[0][j]
                    + self_row[1] * other.m[1][j]
                    + self_row[2] * other.m[2][j];
            }
        }
        Self { m }
    }
}

impl Mat3 {
    /// Construct from rows.
    pub fn from_rows(rows: [[f64; 3]; 3]) -> Self {
        Self { m: rows }
    }

    /// Identity matrix.
    pub fn identity() -> Self {
        Self {
            m: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        }
    }

    /// Right-handed active rotation about the x-axis by `angle` radians.
    pub fn rot_x(angle: f64) -> Self {
        let (s, c) = angle.sin_cos();
        Self {
            m: [[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]],
        }
    }

    /// Right-handed active rotation about the y-axis by `angle` radians.
    pub fn rot_y(angle: f64) -> Self {
        let (s, c) = angle.sin_cos();
        Self {
            m: [[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]],
        }
    }

    /// Right-handed active rotation about the z-axis by `angle` radians.
    pub fn rot_z(angle: f64) -> Self {
        let (s, c) = angle.sin_cos();
        Self {
            m: [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]],
        }
    }

    /// Transpose.
    pub fn transpose(self) -> Self {
        let m = self.m;
        Self {
            m: [
                [m[0][0], m[1][0], m[2][0]],
                [m[0][1], m[1][1], m[2][1]],
                [m[0][2], m[1][2], m[2][2]],
            ],
        }
    }

    /// Apply the matrix to a column vector.
    pub fn apply(self, v: Vec3) -> Vec3 {
        Vec3 {
            xyz: [
                self.m[0][0] * v.xyz[0] + self.m[0][1] * v.xyz[1] + self.m[0][2] * v.xyz[2],
                self.m[1][0] * v.xyz[0] + self.m[1][1] * v.xyz[1] + self.m[1][2] * v.xyz[2],
                self.m[2][0] * v.xyz[0] + self.m[2][1] * v.xyz[1] + self.m[2][2] * v.xyz[2],
            ],
        }
    }

    /// Determinant, expanded along the first row as written.
    pub fn determinant(self) -> f64 {
        let m = self.m;
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    }

    /// Inverse via the adjugate matrix; singular inputs are an error.
    pub fn inverse(self) -> Result<Self> {
        let det = self.determinant();
        if det == 0.0 {
            return Err(Error::SingularMatrix);
        }
        let m = self.m;
        let adj = [
            [
                m[1][1] * m[2][2] - m[1][2] * m[2][1],
                m[0][2] * m[2][1] - m[0][1] * m[2][2],
                m[0][1] * m[1][2] - m[0][2] * m[1][1],
            ],
            [
                m[1][2] * m[2][0] - m[1][0] * m[2][2],
                m[0][0] * m[2][2] - m[0][2] * m[2][0],
                m[0][2] * m[1][0] - m[0][0] * m[1][2],
            ],
            [
                m[1][0] * m[2][1] - m[1][1] * m[2][0],
                m[0][1] * m[2][0] - m[0][0] * m[2][1],
                m[0][0] * m[1][1] - m[0][1] * m[1][0],
            ],
        ];
        let mut inv = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                inv[i][j] = adj[i][j] / det;
            }
        }
        Ok(Self { m: inv })
    }

    /// Maximum absolute deviation of `self^T * self` from the identity;
    /// zero (up to rounding) for exact rotations.
    pub fn orthogonality_error(self) -> f64 {
        let product = self.transpose() * self;
        product
            .m
            .iter()
            .enumerate()
            .flat_map(|(i, row)| {
                row.iter().enumerate().map(move |(j, &v)| {
                    let want = if i == j { 1.0 } else { 0.0 };
                    (v - want).abs()
                })
            })
            .fold(0.0_f64, f64::max)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    const TOL: f64 = 1e-12;

    #[test]
    fn compensated_norm_preserves_residues_and_extreme_domains() {
        let x = f64::from_bits(0x3fee30d883b6a12f);
        let y = f64::from_bits(0x3fd376e83e5e5e19);
        // Independently source-acquired norm, also checked with a 180-digit
        // decimal square/sum/root in the owning Besselian regression.
        let expected = f64::from_bits(0x3fefb87e11fb9479);
        for exponent in [-1000, -40, 0, 1000] {
            let scale = 2.0_f64.powi(exponent);
            for (a, b) in [(x, y), (y, x), (-x, y), (x, -y)] {
                assert_eq!(
                    compensated_hypot(a * scale, b * scale).to_bits(),
                    (expected * scale).to_bits()
                );
            }
        }
        assert_eq!(compensated_hypot(-0.0, 0.0).to_bits(), 0);
        assert_eq!(
            compensated_hypot(f64::from_bits(3), f64::from_bits(4)).to_bits(),
            5
        );
        assert_eq!(compensated_hypot(f64::MAX, 0.0), f64::MAX);
        assert!(compensated_hypot(f64::MAX, f64::MAX).is_infinite());
        assert_eq!(compensated_hypot(f64::NAN, f64::INFINITY), f64::INFINITY);
        assert!(compensated_hypot(f64::NAN, 1.0).is_nan());
    }

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= TOL * b.abs().max(1.0)
    }

    fn close_vec(a: Vec3, b: Vec3) -> bool {
        close(a.xyz[0], b.xyz[0]) && close(a.xyz[1], b.xyz[1]) && close(a.xyz[2], b.xyz[2])
    }

    #[test]
    fn dot_cross_non_axis_aligned() {
        // Independently derived: (1,2,3).(4,5,6) = 32;
        // (1,2,3)x(4,5,6) = (-3,6,-3).
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(4.0, 5.0, 6.0);
        assert!(close(a.dot(b), 32.0));
        assert!(close_vec(a.cross(b), Vec3::new(-3.0, 6.0, -3.0)));
        // Antisymmetry and self-cross.
        assert!(close_vec(b.cross(a), Vec3::new(3.0, -6.0, 3.0)));
        assert_eq!(a.cross(a), Vec3::zero());
    }

    #[test]
    fn norm_and_normalize() {
        let v = Vec3::new(1.0, 2.0, 2.0);
        assert!(close(v.norm(), 3.0));
        let u = v.normalized().unwrap();
        assert!(close(u.norm(), 1.0));
        assert!(close_vec(u, Vec3::new(1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0)));
        assert!(matches!(
            Vec3::zero().normalized(),
            Err(Error::ZeroVectorNormalize)
        ));
    }

    #[test]
    fn rotation_conventions() {
        // Active right-handed rotations: x-axis stays fixed under rot_x,
        // rot_z(90deg) sends +x to +y, rot_y(90deg) sends +z to +x.
        let x = Vec3::new(1.0, 0.0, 0.0);
        let z = Vec3::new(0.0, 0.0, 1.0);
        assert!(close_vec(
            Mat3::rot_x(1.234).apply(x),
            Vec3::new(1.0, 0.0, 0.0)
        ));
        assert!(close_vec(
            Mat3::rot_z(PI / 2.0).apply(x),
            Vec3::new(0.0, 1.0, 0.0)
        ));
        assert!(close_vec(
            Mat3::rot_y(PI / 2.0).apply(z),
            Vec3::new(1.0, 0.0, 0.0)
        ));
        // Non-axis-aligned: rot_z(180deg) negates x and y, keeps z.
        let v = Vec3::new(1.0, 2.0, 3.0);
        assert!(close_vec(
            Mat3::rot_z(PI).apply(v),
            Vec3::new(-1.0, -2.0, 3.0)
        ));
    }

    #[test]
    fn transpose_inverse_and_orthogonality() {
        // det == 1 matrix with independently derived adjugate inverse.
        let a = Mat3::from_rows([[1.0, 2.0, 3.0], [0.0, 1.0, 4.0], [5.0, 6.0, 0.0]]);
        assert!(close(a.determinant(), 1.0));
        let inv = a.inverse().unwrap();
        let want = Mat3::from_rows([[-24.0, 18.0, 5.0], [20.0, -15.0, -4.0], [-5.0, 4.0, 1.0]]);
        for i in 0..3 {
            for j in 0..3 {
                assert!(close(inv.m[i][j], want.m[i][j]), "entry [{i}][{j}]");
            }
        }
        // a * a^-1 is the identity (checks Mul + inverse together).
        let id = a * inv;
        for i in 0..3 {
            for j in 0..3 {
                let want = if i == j { 1.0 } else { 0.0 };
                assert!(close(id.m[i][j], want), "product [{i}][{j}]");
            }
        }
        // Rotations are orthogonal to rounding precision.
        for r in [Mat3::rot_x(0.7), Mat3::rot_y(-1.3), Mat3::rot_z(2.1)] {
            assert!(r.orthogonality_error() < 1e-15);
            // For rotations, transpose is the inverse.
            let v = Vec3::new(1.0, -2.0, 0.5);
            assert!(close_vec(r.transpose().apply(r.apply(v)), v));
        }
    }

    #[test]
    fn singular_matrix_is_an_error() {
        let s = Mat3::from_rows([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]);
        assert_eq!(s.determinant(), 0.0);
        assert!(matches!(s.inverse(), Err(Error::SingularMatrix)));
    }
}
