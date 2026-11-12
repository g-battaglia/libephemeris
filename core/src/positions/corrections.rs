//! Corrections implementation for the parent module.
//! Public entry points and types are exposed through the parent API.

use super::*;

type FrameMatrix = [[f64; 3]; 3];

/// Convert a TT Julian Day to Julian centuries past J2000.
pub(super) fn julian_centuries(jd_tt: f64) -> f64 {
    (jd_tt - FRAME_J2000) / 36525.0
}

/// IAU 2006 mean obliquity in degrees at a finite TT Julian Day.
/// The polynomial uses TT Julian centuries from J2000. Non-finite input is
/// [`Error::NonFiniteArgument`].
pub fn mean_obliquity_iau2006_deg(jd_tt: f64) -> Result<f64> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let t = julian_centuries(jd_tt);
    let mut eps_arcsec = MEAN_OBLIQUITY_COEFFS[0];
    let mut t_power = t;
    for coeff in MEAN_OBLIQUITY_COEFFS.iter().skip(1) {
        eps_arcsec += coeff * t_power;
        t_power *= t;
    }
    Ok(eps_arcsec / 3600.0)
}

/// IAU 2006 general precession in longitude from J2000, in degrees.
/// The argument is a finite TT Julian Day; non-finite input is an error.
pub fn general_precession_deg(jd_tt: f64) -> Result<f64> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let centuries = julian_centuries(jd_tt);
    let mut power = centuries;
    let mut arcseconds = 0.0;
    for coefficient in GENERAL_PRECESSION_COEFFS {
        arcseconds += coefficient * power;
        power *= centuries;
    }
    Ok(arcseconds / 3600.0)
}

/// Derivative of the IAU 2006 longitude-precession polynomial, in degrees
/// per TT day. The argument is a finite TT Julian Day.
pub fn general_precession_rate_deg_day(jd_tt: f64) -> Result<f64> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let t = julian_centuries(jd_tt);
    let rate_arcsec_cy = GENERAL_PRECESSION_COEFFS[0] + 2.0 * GENERAL_PRECESSION_COEFFS[1] * t;
    Ok(rate_arcsec_cy / (3600.0 * 36525.0))
}

/// IAU 2006 Fukushima-Williams precession angles (IERS 2010, Table 5.1).
///
/// Returns `(gamb, phib, psib, epsa)` in radians at a finite TT Julian Day.
/// Non-finite input is [`Error::NonFiniteArgument`].
pub fn iau2006_precession_angles(jd_tt: f64) -> Result<(f64, f64, f64, f64)> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let t = julian_centuries(jd_tt);
    let gamb = (-0.052928
        + (10.556378
            + (0.4932044 + (-0.00031238 + (-0.000002788 + 0.0000000260 * t) * t) * t) * t)
            * t)
        * ARCSEC_TO_RAD;
    let phib = (84381.412819
        + (-46.811016
            + (0.0511268 + (0.00053289 + (-0.000000440 - 0.0000000176 * t) * t) * t) * t)
            * t)
        * ARCSEC_TO_RAD;
    let psib = (-0.041775
        + (5038.481484
            + (1.5584175 + (-0.00018522 + (-0.000026452 - 0.0000000148 * t) * t) * t) * t)
            * t)
        * ARCSEC_TO_RAD;
    let epsa = (84381.406
        + (-46.836769
            + (-0.0001831 + (0.00200340 + (-0.000000576 - 0.0000000434 * t) * t) * t) * t)
            * t)
        * ARCSEC_TO_RAD;
    Ok((gamb, phib, psib, epsa))
}

/// Fukushima-Williams bias/precession/nutation rotation from radian angles.
///
/// Return a row-major matrix acting on a Cartesian column vector:
/// `R = R1(-eps) R3(-psi) R1(phib) R3(gamb)`. `gamb` and `phib` are
/// the frame-bias/precession orientation angles; `psi` is precession in
/// longitude (plus nutation when requested), and `eps` the corresponding
/// mean or true obliquity. With IAU 2006/2000A angles, the map is GCRS to
/// the selected of-date equatorial frame. Epoch/model validity belongs to
/// the supplying angle routine; this function has no time argument.
/// Execute ordered pre-rotations with paired trigonometry and fused final
/// row updates, avoiding regrouped products that change small CIP entries.
/// Non-finite inputs or resulting matrix entries are refused.
pub fn fw2m(gamb: f64, phib: f64, psi: f64, eps: f64) -> Result<[[f64; 3]; 3]> {
    for v in [gamb, phib, psi, eps] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    let matrix = crate::erfa::fw2m_rotation(gamb, phib, psi, eps);
    for row in matrix {
        if let Some(&value) = row.iter().find(|value| !value.is_finite()) {
            return Err(Error::MathDomain { value });
        }
    }
    Ok(matrix)
}

/// Multiply a dimensionless row-major 3x3 matrix by a finite vector.
/// The output has the vector's units. Nonfinite inputs return a domain error.
pub fn mat3_vec3(m: [[f64; 3]; 3], v: [f64; 3]) -> Result<[f64; 3]> {
    check_finite_vec3(v)?;
    for row in m {
        check_finite_vec3(row)?;
    }
    let result = m.map(|row| compensated_dot3(row, v));
    check_finite_vec3(result)?;
    Ok(result)
}

/// Three-product dot with compensated product and addition residuals.
/// Fused multiply-add recovers the rounding error of each product; Neumaier
/// compensation retains small terms in cancellation-prone rotation matrices.
/// All operations use binary64. This kernel is shared by matrix and vector
/// products and is independent of platform linear-algebra libraries.
pub(crate) fn compensated_dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    let mut sum = 0.0_f64;
    let mut correction = 0.0_f64;
    for (left, right) in a.into_iter().zip(b) {
        let product = left * right;
        let residual = left.mul_add(right, -product);
        let next = sum + product;
        let addition_error = if sum.abs() >= product.abs() {
            (sum - next) + product
        } else {
            (product - next) + sum
        };
        correction += residual + addition_error;
        sum = next;
    }
    sum + correction
}

/// Matrix product `C = A @ B`, inner accumulation in `k` order.
pub fn mat3_mul(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> Result<[[f64; 3]; 3]> {
    for row in a {
        check_finite_vec3(row)?;
    }
    for row in b {
        check_finite_vec3(row)?;
    }
    let mut m = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            m[i][j] = compensated_dot3(a[i], [b[0][j], b[1][j], b[2][j]]);
            if !m[i][j].is_finite() {
                return Err(Error::MathDomain { value: m[i][j] });
            }
        }
    }
    Ok(m)
}

/// Rotate an equatorial Cartesian vector into the ecliptic frame.
///
/// `eps_rad` is the obliquity in radians; Cartesian units are unchanged.
/// The common compensated matrix product limits cancellation near the
/// ecliptic plane. Nonfinite inputs and arithmetic overflow are rejected.
pub fn rotate_equatorial_to_ecliptic(v: [f64; 3], eps_rad: f64) -> Result<[f64; 3]> {
    check_finite_vec3(v)?;
    if !eps_rad.is_finite() {
        return Err(Error::NonFiniteArgument { value: eps_rad });
    }
    let (sin_eps, cos_eps) = eps_rad.sin_cos();
    mat3_vec3(
        [
            [1.0, 0.0, 0.0],
            [0.0, cos_eps, sin_eps],
            [0.0, -sin_eps, cos_eps],
        ],
        v,
    )
}

/// ICRS to IAU 2006 mean-ecliptic J2000 rotation, including frame bias.
/// The dimensionless matrix is `Rx(obliquity) * frame_bias`.
pub fn icrs_to_ecl_j2000_matrix() -> [[f64; 3]; 3] {
    let ce = OBLIQUITY_J2000_REF_RAD.cos();
    let se = OBLIQUITY_J2000_REF_RAD.sin();
    let rx = [[1.0, 0.0, 0.0], [0.0, ce, se], [0.0, -se, ce]];
    std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            compensated_dot3(rx[i], std::array::from_fn(|k| FRAME_BIAS_J2000[k][j]))
        })
    })
}

/// Rotate a finite ICRS vector to the mean equator/equinox of J2000.
/// The frame-bias rotation preserves the input vector's length unit.
pub fn rotate_icrs_to_j2000_mean_equator(v: [f64; 3]) -> Result<[f64; 3]> {
    mat3_vec3(FRAME_BIAS_J2000, v)
}

/// Rotate a finite ICRS vector to the IAU 2006 mean ecliptic of J2000.
/// The rotation preserves the input vector's length unit.
pub fn rotate_icrs_to_ecliptic_j2000(v: [f64; 3]) -> Result<[f64; 3]> {
    mat3_vec3(icrs_to_ecl_j2000_matrix(), v)
}

/// Convert a finite Cartesian vector to longitude, latitude and radius.
///
/// Angles are degrees; longitude is normalized to `[0, 360)`. Radius retains
/// the input length unit. The zero vector returns three zeros.
pub fn cartesian_to_spherical_deg(v: [f64; 3]) -> Result<[f64; 3]> {
    check_finite_vec3(v)?;
    let r = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if r == 0.0 {
        return Ok([0.0, 0.0, 0.0]);
    }
    let mut lon = v[1].atan2(v[0]).to_degrees() % 360.0;
    if lon < 0.0 {
        lon += 360.0;
    }
    let lat = (v[2] / r).clamp(-1.0, 1.0).asin().to_degrees();
    Ok([lon, lat, r])
}

/// Rotate a direction between ecliptic and equatorial coordinates.
///
/// Longitude, latitude and obliquity `eps_deg` are degrees. Negative obliquity
/// maps ecliptic longitude/latitude to right ascension/declination; positive
/// obliquity maps back. Both entries apply an x-axis rotation by `-eps_deg`
/// to the unit direction, then recover longitude with `atan2(y, x)` and
/// latitude with `asin(z / r)`, where `r` is the rotated vector's length.
/// Longitude is reduced modulo 360 degrees; zero has the positive sign, but
/// floating-point rounding can produce 360 degrees at the upper endpoint.
/// Latitude retains its computed sign. Finite angles outside one turn or
/// the conventional latitude interval are evaluated through their full
/// trigonometric components. Any non-finite argument returns
/// [`Error::NonFiniteArgument`]. This closed-form operation performs no I/O.
pub fn cotrans(lon_deg: f64, lat_deg: f64, eps_deg: f64) -> Result<(f64, f64)> {
    for v in [lon_deg, lat_deg, eps_deg] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    let eps_rad = (-eps_deg).to_radians();
    let lon_rad = lon_deg.to_radians();
    let lat_rad = lat_deg.to_radians();
    let cos_eps = eps_rad.cos();
    let sin_eps = eps_rad.sin();
    let cos_lat = lat_rad.cos();
    let sin_lat = lat_rad.sin();
    let cos_lon = lon_rad.cos();
    let sin_lon = lon_rad.sin();
    let x = cos_lat * cos_lon;
    let y = cos_lat * sin_lon * cos_eps - sin_lat * sin_eps;
    let z = cos_lat * sin_lon * sin_eps + sin_lat * cos_eps;
    let mut new_lon = y.atan2(x).to_degrees() % 360.0;
    if new_lon < 0.0 {
        new_lon += 360.0;
    }
    // Canonicalize only longitude zero, preserving nonzero remainders.
    if new_lon == 0.0 {
        new_lon = 0.0;
    }
    let r = (x * x + y * y + z * z).sqrt();
    let new_lat = if r > 0.0 {
        (z / r).clamp(-1.0, 1.0).asin().to_degrees()
    } else {
        0.0
    };
    Ok((new_lon, new_lat))
}

/// Mean-to-true equatorial nutation rotation matrix:
/// `N = R1(-(eps + deps)) @ R3(-dpsi) @ R1(eps)`.
///
/// `eps_mean_rad` is mean obliquity, `dpsi` nutation in longitude and
/// `deps` nutation in obliquity, all in radians. The row-major matrix
/// acts on an equatorial column vector without changing its units.
/// The native ERFA rotation sequence applies the rightmost rotation
/// first, retaining paired trigonometry and fused row updates rather
/// than regrouping matrix products. Zero nutation gives the identity
/// up to rounding. Non-finite inputs or matrix entries are refused.
pub fn nutation_matrix(eps_mean_rad: f64, dpsi: f64, deps: f64) -> Result<[[f64; 3]; 3]> {
    for v in [eps_mean_rad, dpsi, deps] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    let n = crate::erfa::numat_rotation(eps_mean_rad, dpsi, deps);
    for row in n {
        check_finite_vec3(row)?;
    }
    Ok(n)
}

/// Vondrák et al. (2011) mean-ecliptic pole at the TT Julian Day.
pub(crate) fn ecliptic_pole_vondrak(jd_tt: f64) -> [f64; 3] {
    ecliptic_pole_vondrak_from_centuries(julian_centuries(jd_tt))
}

/// Shared pole series over its already acquired direct century operand.
fn ecliptic_pole_vondrak_from_centuries(t: f64) -> [f64; 3] {
    let two_pi = 2.0 * std::f64::consts::PI;
    let mut p = 0.0;
    let mut q = 0.0;
    for i in 0..VONDRAK_PQPER_PERIODS.len() {
        let a = two_pi * t / VONDRAK_PQPER_PERIODS[i];
        let (s, c) = a.sin_cos();
        p += c * VONDRAK_PQPER_P_COS[i] + s * VONDRAK_PQPER_P_SIN[i];
        q += c * VONDRAK_PQPER_Q_COS[i] + s * VONDRAK_PQPER_Q_SIN[i];
    }
    let mut w = 1.0;
    for poly in VONDRAK_PQPOL {
        p += poly.0 * w;
        q += poly.1 * w;
        w *= t;
    }
    p *= ARCSEC_TO_RAD;
    q *= ARCSEC_TO_RAD;
    let mut z = 1.0 - p * p - q * q;
    z = if z > 0.0 { z.sqrt() } else { 0.0 };
    let se = VONDRAK_EPS0_RAD.sin();
    let ce = VONDRAK_EPS0_RAD.cos();
    [p, -q * ce - z * se, -q * se + z * ce]
}

/// Vondrák et al. (2011) mean-equator pole at the TT Julian Day.
pub(crate) fn equator_pole_vondrak(jd_tt: f64) -> [f64; 3] {
    equator_pole_vondrak_from_centuries(julian_centuries(jd_tt))
}

/// Vondrák mean-equator pole series as a function of TT centuries from J2000.
fn equator_pole_vondrak_from_centuries(t: f64) -> [f64; 3] {
    let two_pi = 2.0 * std::f64::consts::PI;
    let mut x = 0.0;
    let mut y = 0.0;
    for i in 0..VONDRAK_XYPER_PERIODS.len() {
        let a = two_pi * t / VONDRAK_XYPER_PERIODS[i];
        let (s, c) = a.sin_cos();
        x += c * VONDRAK_XYPER_X_COS[i] + s * VONDRAK_XYPER_X_SIN[i];
        y += c * VONDRAK_XYPER_Y_COS[i] + s * VONDRAK_XYPER_Y_SIN[i];
    }
    let mut w = 1.0;
    for poly in VONDRAK_XYPOL {
        x += poly.0 * w;
        y += poly.1 * w;
        w *= t;
    }
    x *= ARCSEC_TO_RAD;
    y *= ARCSEC_TO_RAD;
    let w2 = x * x + y * y;
    let z = if w2 < 1.0 { (1.0 - w2).sqrt() } else { 0.0 };
    [x, y, z]
}

/// Evaluate the direct Vondrák et al. (2011) mean-obliquity series in radians.
///
/// This diagnostic series is distinct from the pole-angle obliquity used to
/// keep the native precession and equator/ecliptic transformations consistent.
/// The argument is a finite TT Julian Day.
pub fn vondrak_mean_obliquity_series_rad(jd_tt: f64) -> Result<f64> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let t = julian_centuries(jd_tt);
    let two_pi = 2.0 * std::f64::consts::PI;
    let mut q = 0.0;
    for i in 0..VONDRAK_PEPER_PERIODS.len() {
        let a = two_pi * t / VONDRAK_PEPER_PERIODS[i];
        q += a.cos() * VONDRAK_PEPER_Q_COS[i] + a.sin() * VONDRAK_PEPER_Q_SIN[i];
    }
    let mut w = 1.0;
    for poly in VONDRAK_PEPOL {
        q += poly.1 * w;
        w *= t;
    }
    Ok(q * ARCSEC_TO_RAD)
}

/// Mean obliquity in radians from the angle between the Vondrák et al.
/// (2011) mean-equator and mean-ecliptic poles at a finite TT Julian Day.
pub fn vondrak_mean_obliquity_rad(jd_tt: f64) -> Result<f64> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let e = ecliptic_pole_vondrak(jd_tt);
    let q = equator_pole_vondrak(jd_tt);
    // Bound the direction cosine to the domain of acos despite roundoff.
    let d = (e[0] * q[0] + e[1] * q[1] + e[2] * q[2]).clamp(-1.0, 1.0);
    Ok(d.acos())
}

/// Mean obliquity in radians from Vondrák poles at TT centuries from J2000.
///
/// This low-level entry evaluates `(JD_TT − 2451545) / 36525` directly, without
/// an additional epoch conversion. Non-finite arguments propagate through
/// ordinary binary64 arithmetic. The published long-term model's physical
/// uncertainty grows away from J2000; finite output is not an accuracy bound.
pub fn source_vondrak_mean_obliquity_from_centuries(t: f64) -> f64 {
    let e = ecliptic_pole_vondrak_from_centuries(t);
    let q = equator_pole_vondrak_from_centuries(t);
    let d = (e[0] * q[0] + e[1] * q[1] + e[2] * q[2]).clamp(-1.0, 1.0);
    d.acos()
}

/// Mean obliquity in degrees at a finite TT Julian Day, from Vondrák poles.
pub fn vondrak_mean_obliquity_deg(jd_tt: f64) -> Result<f64> {
    Ok(vondrak_mean_obliquity_rad(jd_tt)?.to_degrees())
}

/// Vondrák 2011/2012 J2000-mean to of-date mean-equatorial precession.
///
/// `jd_tt` is a TT Julian Day. The row-major matrix maps a J2000
/// mean-equatorial column vector to the mean equator/equinox of date,
/// with rows `(equinox, middle, equator pole)` and no ICRS frame bias.
/// Convert first to Julian epoch `2000 + (jd_tt - J2000) / 365.25`,
/// then evaluate the native ERFA pole series and ordered cross-product
/// normalization. This acquisition is distinct from the scalar-JD pole
/// series used for mean obliquity. The published model spans ±200,000
/// years about J2000, with growing physical uncertainty, not certified
/// astronomical accuracy. Non-finite epochs, arithmetic overflow and
/// an undefined equinox are refused.
pub fn vondrak_precession_matrix_nobias(jd_tt: f64) -> Result<[[f64; 3]; 3]> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    crate::erfa::ltp_matrix(crate::erfa::julian_epoch(jd_tt)?)
}

/// Mean-equator-of-date precession, optionally including ICRS frame bias.
///
/// `jd_tt` is a finite TT Julian Day. The row-major matrix acts on Cartesian
/// column vectors. With bias, the IERS 2010 small-angle realization maps ICRS
/// to date; without it, the input is the J2000 mean frame. Nutation is excluded.
/// Domains, errors and model uncertainty follow
/// [`vondrak_precession_matrix_nobias`].
pub fn vondrak_precession_matrix(jd_tt: f64, frame_bias: bool) -> Result<[[f64; 3]; 3]> {
    let epj = crate::erfa::julian_epoch(jd_tt)?;
    if frame_bias {
        crate::erfa::ltpb_matrix(epj)
    } else {
        crate::erfa::ltp_matrix(epj)
    }
}

/// ICRS or J2000-mean input to true-equatorial-of-date rotation.
///
/// `jd_tt` is TT JD; nutation angles are radians. Return `(N @ P, eps_true)`
/// with row-major matrices and `eps_true = eps_mean + deps` in radians.
/// `frame_bias` selects ICRS rather than J2000 mean input. Mean obliquity uses
/// the Vondrák pole angle; precession uses its Julian-epoch realization.
/// Non-finite inputs and component frame errors propagate. Long-term model
/// uncertainty remains that of the published precession series.
pub fn vondrak_pn_matrix(
    jd_tt: f64,
    dpsi_rad: f64,
    deps_rad: f64,
    frame_bias: bool,
) -> Result<([[f64; 3]; 3], f64)> {
    let eps_mean = vondrak_mean_obliquity_rad(jd_tt)?;
    for v in [dpsi_rad, deps_rad] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    let n = nutation_matrix(eps_mean, dpsi_rad, deps_rad)?;
    let p = vondrak_precession_matrix(jd_tt, frame_bias)?;
    Ok((pn_matrix_product(n, p)?, eps_mean + deps_rad))
}

/// Compose dimensionless nutation and precession rotations.
/// Nonfinite inputs are rejected by the common matrix kernel.
pub(super) fn pn_matrix_product(n: [[f64; 3]; 3], p: [[f64; 3]; 3]) -> Result<[[f64; 3]; 3]> {
    mat3_mul(n, p)
}

impl FrameMatrices {
    /// Checked constructor: every entry and angle must be finite.
    pub fn new(
        p: [[f64; 3]; 3],
        p_nobias: [[f64; 3]; 3],
        pn: [[f64; 3]; 3],
        pn_nobias: [[f64; 3]; 3],
        eps_true_rad: f64,
        eps_mean_rad: f64,
    ) -> Result<Self> {
        for m in [p, p_nobias, pn, pn_nobias] {
            for row in m {
                check_finite_vec3(row)?;
            }
        }
        for v in [eps_true_rad, eps_mean_rad] {
            if !v.is_finite() {
                return Err(Error::NonFiniteArgument { value: v });
            }
        }
        Ok(Self {
            p,
            p_nobias,
            pn,
            pn_nobias,
            eps_true_rad,
            eps_mean_rad,
            native_ecliptic: None,
        })
    }

    /// Biased ICRS -> mean-equator-of-date matrix.
    pub fn p(self) -> [[f64; 3]; 3] {
        self.p
    }

    /// Bias-free mean-equator-of-date matrix (`FLG_ICRS` companion).
    pub fn p_nobias(self) -> [[f64; 3]; 3] {
        self.p_nobias
    }

    /// Biased ICRS -> true-equator-of-date matrix.
    pub fn pn(self) -> [[f64; 3]; 3] {
        self.pn
    }

    /// Bias-free true-equator-of-date matrix (`FLG_ICRS` companion).
    pub fn pn_nobias(self) -> [[f64; 3]; 3] {
        self.pn_nobias
    }

    /// True obliquity of date in radians (bias-independent).
    pub fn eps_true_rad(self) -> f64 {
        self.eps_true_rad
    }

    /// Mean obliquity of date in radians (Vondrák pole-angle value).
    pub fn eps_mean_rad(self) -> f64 {
        self.eps_mean_rad
    }
}

/// Native Vondrák frame bundle for `jd_tt` with caller-supplied nutation
/// angles (radians; LEB Chebyshev nutation on the fast path).
pub fn frame_matrices_for_date(jd_tt: f64, dpsi_rad: f64, deps_rad: f64) -> Result<FrameMatrices> {
    let key = [jd_tt.to_bits(), dpsi_rad.to_bits(), deps_rad.to_bits()];
    if let Some(frames) = crate::pure_cache::frame(key) {
        return Ok(frames);
    }
    let frames = frame_matrices_uncached(jd_tt, dpsi_rad, deps_rad)?;
    crate::pure_cache::retain_frame(key, frames);
    Ok(frames)
}

fn frame_matrices_uncached(jd_tt: f64, dpsi_rad: f64, deps_rad: f64) -> Result<FrameMatrices> {
    // Nutation and obliquity depend on the epoch, not on frame bias. Share
    // them without changing the ordered products for either precession frame.
    let eps_mean = vondrak_mean_obliquity_rad(jd_tt)?;
    for value in [dpsi_rad, deps_rad] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    let nutation = nutation_matrix(eps_mean, dpsi_rad, deps_rad)?;
    let p = vondrak_precession_matrix(jd_tt, true)?;
    let p_nobias = vondrak_precession_matrix(jd_tt, false)?;
    let pn = pn_matrix_product(nutation, p)?;
    let pn_nobias = pn_matrix_product(nutation, p_nobias)?;
    let eps_true = eps_mean + deps_rad;
    let mut frames = FrameMatrices::new(p, p_nobias, pn, pn_nobias, eps_true, eps_mean)?;
    frames.native_ecliptic = Some((mean_ecliptic_basis(jd_tt)?, dpsi_rad));
    Ok(frames)
}

/// Mean-ecliptic orthonormal basis from the published Vondrák TT poles.
/// The equinox is the normalized cross product of equator and ecliptic poles;
/// the second row completes the right-handed basis. All rows are dimensionless.
/// Direct TT centuries avoid losing fractional-day precision through a large
/// Julian-year epoch. Nonfinite poles or an undefined equinox are refused.
fn mean_ecliptic_basis(jd_tt: f64) -> Result<[[f64; 3]; 3]> {
    let ecliptic = ecliptic_pole_vondrak(jd_tt);
    let equator = equator_pole_vondrak(jd_tt);
    check_finite_vec3(ecliptic)?;
    check_finite_vec3(equator)?;
    let cross = |a: [f64; 3], b: [f64; 3]| {
        std::array::from_fn(|i| {
            let j = (i + 1) % 3;
            let k = (i + 2) % 3;
            compensated_dot3([a[j], -a[k], 0.0], [b[k], b[j], 0.0])
        })
    };
    let equinox = cross(equator, ecliptic);
    let length = vec3_dist(equinox);
    if length == 0.0 || !length.is_finite() {
        return Err(Error::MathDomain { value: length });
    }
    let equinox = equinox.map(|value| value / length);
    Ok([equinox, cross(ecliptic, equinox), ecliptic])
}

/// The same direct-TT ecliptic basis and its analytic derivative per TT day.
/// Differentiate every retained periodic/polynomial pole term, then the cross
/// product and its normalization. The central basis keeps its scalar arithmetic.
fn vondrak_pole_rates(t: f64, equator: [f64; 3], ecliptic: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let periodic_rate = |periods: &[f64], cosine: &[f64], sine: &[f64]| {
        let mut rate = 0.0;
        for ((&period, &c), &s) in periods.iter().zip(cosine).zip(sine) {
            let frequency = std::f64::consts::TAU / period;
            let (sin, cos) = (std::f64::consts::TAU * t / period).sin_cos();
            rate += frequency * (-sin * c + cos * s);
        }
        rate
    };
    let polynomial_rate = |rows: &[(f64, f64)], component: usize| {
        let mut rate = 0.0;
        let mut power = 1.0;
        for (order, &(a, b)) in rows.iter().enumerate().skip(1) {
            rate += order as f64 * if component == 0 { a } else { b } * power;
            power *= t;
        }
        rate
    };
    let scale = ARCSEC_TO_RAD / 36525.0;
    let pd = (periodic_rate(
        &VONDRAK_PQPER_PERIODS,
        &VONDRAK_PQPER_P_COS,
        &VONDRAK_PQPER_P_SIN,
    ) + polynomial_rate(&VONDRAK_PQPOL, 0))
        * scale;
    let qd = (periodic_rate(
        &VONDRAK_PQPER_PERIODS,
        &VONDRAK_PQPER_Q_COS,
        &VONDRAK_PQPER_Q_SIN,
    ) + polynomial_rate(&VONDRAK_PQPOL, 1))
        * scale;
    let xd = (periodic_rate(
        &VONDRAK_XYPER_PERIODS,
        &VONDRAK_XYPER_X_COS,
        &VONDRAK_XYPER_X_SIN,
    ) + polynomial_rate(&VONDRAK_XYPOL, 0))
        * scale;
    let yd = (periodic_rate(
        &VONDRAK_XYPER_PERIODS,
        &VONDRAK_XYPER_Y_COS,
        &VONDRAK_XYPER_Y_SIN,
    ) + polynomial_rate(&VONDRAK_XYPOL, 1))
        * scale;
    let equator_zd = if equator[2] > 0.0 {
        -(equator[0] * xd + equator[1] * yd) / equator[2]
    } else {
        0.0
    };
    let equator_rate = [xd, yd, equator_zd];
    let (s, c) = VONDRAK_EPS0_RAD.sin_cos();
    let q = -ecliptic[1] * c - ecliptic[2] * s;
    let z = -ecliptic[1] * s + ecliptic[2] * c;
    let zd = if z > 0.0 {
        -(ecliptic[0] * pd + q * qd) / z
    } else {
        0.0
    };
    let ecliptic_rate = [pd, -qd * c - zd * s, -qd * s + zd * c];
    (equator_rate, ecliptic_rate)
}

/// Differentiate the shared equinox normalization and the chosen north pole.
#[allow(clippy::too_many_arguments)]
fn pole_basis_rate(
    basis: [[f64; 3]; 3],
    equator: [f64; 3],
    ecliptic: [f64; 3],
    equator_rate: [f64; 3],
    ecliptic_rate: [f64; 3],
    pole_rate: [f64; 3],
) -> Result<[[f64; 3]; 3]> {
    let cross = |a: [f64; 3], b: [f64; 3]| -> [f64; 3] {
        std::array::from_fn(|i| {
            let j = (i + 1) % 3;
            let k = (i + 2) % 3;
            a[j] * b[k] - a[k] * b[j]
        })
    };
    let add = |a: [f64; 3], b: [f64; 3]| std::array::from_fn(|i| a[i] + b[i]);
    let equinox_rate = add(cross(equator_rate, ecliptic), cross(equator, ecliptic_rate));
    let length = vec3_dist(cross(equator, ecliptic));
    let longitudinal = compensated_dot3(basis[0], equinox_rate);
    let equinox_rate =
        std::array::from_fn(|i| (equinox_rate[i] - basis[0][i] * longitudinal) / length);
    let middle_rate = add(cross(pole_rate, basis[0]), cross(basis[2], equinox_rate));
    let rate = [equinox_rate, middle_rate, pole_rate];
    for row in rate {
        check_finite_vec3(row)?;
    }
    Ok(rate)
}

/// Direct-TT ecliptic basis and analytic derivative in inverse TT days.
pub(super) fn mean_ecliptic_basis_and_rate(jd_tt: f64) -> Result<(FrameMatrix, FrameMatrix)> {
    let basis = mean_ecliptic_basis(jd_tt)?;
    let equator = equator_pole_vondrak(jd_tt);
    let (equator_rate, ecliptic_rate) =
        vondrak_pole_rates(julian_centuries(jd_tt), equator, basis[2]);
    let rate = pole_basis_rate(
        basis,
        equator,
        basis[2],
        equator_rate,
        ecliptic_rate,
        ecliptic_rate,
    )?;
    Ok((basis, rate))
}

/// Precession and its analytic TT rate on the scalar frame's Julian-epoch model.
/// The epoch realization stays distinct from direct-TT obliquity/ecliptic poles.
pub(super) fn precession_matrix_and_rate(
    jd_tt: f64,
    bias: bool,
) -> Result<(FrameMatrix, FrameMatrix)> {
    let epj = crate::erfa::julian_epoch(jd_tt)?;
    let equator = crate::erfa::ltpequ(epj)?;
    let ecliptic = crate::erfa::ltpecl(epj)?;
    let (equator_rate, ecliptic_rate) =
        vondrak_pole_rates((epj - 2000.0) / 100.0, equator, ecliptic);
    let basis = vondrak_precession_matrix(jd_tt, false)?;
    let mut rate = pole_basis_rate(
        basis,
        equator,
        ecliptic,
        equator_rate,
        ecliptic_rate,
        equator_rate,
    )?;
    let central = if bias {
        rate = mat3_mul(rate, FRAME_BIAS_LTPB)?;
        vondrak_precession_matrix(jd_tt, true)?
    } else {
        basis
    };
    Ok((central, rate))
}

/// Mean-ecliptic scalar-epoch basis and its complete TT derivative.
/// The sidereal scalar law first acquires a Julian epoch. Keep that operand
/// realization distinct from the direct-TT apparent-frame basis, and reverse
/// the two equinox-plane rows to match the sidereal ecliptic convention.
fn sidereal_ecliptic_basis_and_rate(jd_tt: f64) -> Result<(FrameMatrix, FrameMatrix)> {
    let epj = crate::erfa::julian_epoch(jd_tt)?;
    let equator = crate::erfa::ltpequ(epj)?;
    let ecliptic = crate::erfa::ltpecl(epj)?;
    let scalar = crate::erfa::ltp_ecliptic_frame(jd_tt)?;
    let basis = [scalar[0], scalar[1], scalar[2]];
    let positive_basis = [basis[0].map(|v| -v), basis[1].map(|v| -v), basis[2]];
    let (equator_rate, ecliptic_rate) =
        vondrak_pole_rates((epj - 2000.0) / 100.0, equator, ecliptic);
    let positive_rate = pole_basis_rate(
        positive_basis,
        equator,
        ecliptic,
        equator_rate,
        ecliptic_rate,
        ecliptic_rate,
    )?;
    let rate = [
        positive_rate[0].map(|v| -v),
        positive_rate[1].map(|v| -v),
        positive_rate[2],
    ];
    Ok((basis, rate))
}

/// Derivative of the selected epoch-pair mean ayanamsha, degrees per TT day.
///
/// `defining_epoch_tt` and `defining_value_deg` are immutable acquired model
/// operands. Method B projects the moving equinox onto the defining ecliptic;
/// ECL_DATE projects the fixed defining zero direction onto the moving
/// ecliptic. Differentiate each actual atan2 through `(x*y'-y*x')/(x²+y²)`.
/// No general-precession polynomial substitutes for either scalar law.
/// Finite input, pole normalization and undefined projections use typed errors.
pub(super) fn epoch_pair_ayanamsha_rate_deg_day(
    jd_tt: f64,
    defining_value_deg: f64,
    defining_epoch_tt: f64,
    ecl_date: bool,
) -> Result<f64> {
    for value in [jd_tt, defining_value_deg, defining_epoch_tt] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    let fixed = crate::erfa::ltp_ecliptic_frame(defining_epoch_tt)?;
    let (moving, rate) = sidereal_ecliptic_basis_and_rate(jd_tt)?;
    let (x, y, xd, yd, sign) = if ecl_date {
        let (sine, cosine) = defining_value_deg.to_radians().sin_cos();
        let zero = std::array::from_fn(|i| fixed[0][i] * cosine + fixed[1][i] * sine);
        (
            compensated_dot3(moving[0], zero),
            compensated_dot3(moving[1], zero),
            compensated_dot3(rate[0], zero),
            compensated_dot3(rate[1], zero),
            1.0,
        )
    } else {
        (
            compensated_dot3(fixed[0], moving[0]),
            compensated_dot3(fixed[1], moving[0]),
            compensated_dot3(fixed[0], rate[0]),
            compensated_dot3(fixed[1], rate[0]),
            -1.0,
        )
    };
    let radius_squared = x.mul_add(x, y * y);
    if radius_squared == 0.0 || !radius_squared.is_finite() {
        return Err(Error::MathDomain {
            value: radius_squared,
        });
    }
    let numerator = compensated_dot3([x, -y, 0.0], [yd, xd, 0.0]);
    let derivative = (sign * numerator / radius_squared).to_degrees();
    if !derivative.is_finite() {
        return Err(Error::MathDomain { value: derivative });
    }
    Ok(derivative)
}

/// Pole-angle obliquity and its rate, radians and radians per TT day.
pub(super) fn mean_obliquity_and_rate(jd_tt: f64) -> Result<(f64, f64)> {
    let angle = vondrak_mean_obliquity_rad(jd_tt)?;
    let equator = equator_pole_vondrak(jd_tt);
    let ecliptic = ecliptic_pole_vondrak(jd_tt);
    let (qd, ed) = vondrak_pole_rates(julian_centuries(jd_tt), equator, ecliptic);
    let cosine = (ecliptic[0] * equator[0] + ecliptic[1] * equator[1] + ecliptic[2] * equator[2])
        .clamp(-1.0, 1.0);
    let rate = -(compensated_dot3(qd, ecliptic) + compensated_dot3(equator, ed))
        / (1.0 - cosine * cosine).sqrt();
    if !rate.is_finite() {
        return Err(Error::MathDomain { value: rate });
    }
    Ok((angle, rate))
}

/// Rotate an apparent ICRS position in AU into the requested output frame.
///
/// `jd_tt` is TT Julian Day. `frames` supplies biased and bias-free matrices
/// plus obliquities for that epoch. The position and central-difference speed
/// paths use the same transformation, re-evaluated at each sample epoch to
/// include precession/nutation motion; the finite-difference speed is not an
/// exact analytic derivative. `flags` selects the branches in this order:
///
/// - `EQUATORIAL|J2000`: mean equator/equinox of J2000 (bias, unless `ICRS`).
/// - `EQUATORIAL|SIDEREAL`: mean equator of date, no nutation (bias, unless
///   `ICRS`).
/// - `EQUATORIAL`: true equator of date, or mean when `NONUT` (bias, unless
///   `ICRS`, which rebuilds the bias-free true matrix from the same nutation).
/// - `J2000`: mean ecliptic of J2000 (bias, unless `ICRS`, which rotates by
///   the J2000 mean obliquity only).
/// - default: true ecliptic of date, or mean when `NONUT` (`ICRS` drops the
///   bias; the true obliquity is bias-independent).
///
/// Native ecliptic-of-date frames project onto the Vondrák ecliptic basis,
/// then apply longitude nutation. The identity
/// `Rx(eps+deps) N(eps,dpsi,deps) P = Rz(dpsi) Rx(eps) P` makes latitude
/// independent of nutation. This avoids cancellation through a rounded
/// true-equatorial matrix and obliquity. Explicit custom matrices supplied
/// through [`FrameMatrices::new`] retain their own matrix/obliquity route.
/// Returns spherical `(lon_deg, lat_deg, dist_au)` unless `want_xyz`, which
/// returns Cartesian AU in the selected frame without a spherical roundtrip.
/// A zero vector has three zero spherical components. Invalid inputs and
/// frame-model arithmetic failures propagate from the selected components;
/// physical model validity follows the precession/nutation constructors.
pub fn frame_transform(
    geo: [f64; 3],
    jd_tt: f64,
    flags: i32,
    want_xyz: bool,
    frames: FrameMatrices,
) -> Result<[f64; 3]> {
    check_finite_vec3(geo)?;
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let icrs = flags & FLG_ICRS != 0;
    let v = if flags & FLG_EQUATORIAL != 0 && flags & FLG_J2000 != 0 {
        if icrs {
            geo
        } else {
            mat3_vec3(FRAME_BIAS_J2000, geo)?
        }
    } else if flags & FLG_EQUATORIAL != 0 && flags & FLG_SIDEREAL != 0 {
        let m = if icrs { frames.p_nobias } else { frames.p };
        mat3_vec3(m, geo)?
    } else if flags & FLG_EQUATORIAL != 0 {
        if flags & FLG_NONUT != 0 {
            let m = if icrs { frames.p_nobias } else { frames.p };
            mat3_vec3(m, geo)?
        } else if icrs {
            mat3_vec3(frames.pn_nobias, geo)?
        } else {
            mat3_vec3(frames.pn, geo)?
        }
    } else if flags & FLG_J2000 != 0 {
        if icrs {
            let eps_j2000 = vondrak_mean_obliquity_deg(FRAME_J2000)?.to_radians();
            rotate_equatorial_to_ecliptic(geo, eps_j2000)?
        } else {
            mat3_vec3(icrs_to_ecl_j2000_matrix(), geo)?
        }
    } else if let Some((mean_ecliptic, dpsi)) = frames.native_ecliptic {
        let input = if icrs {
            geo
        } else {
            mat3_vec3(FRAME_BIAS_LTPB, geo)?
        };
        let mut ecliptic = mat3_vec3(mean_ecliptic, input)?;
        if flags & FLG_NONUT == 0 {
            let (sine, cosine) = dpsi.sin_cos();
            let [x, y, z] = ecliptic;
            ecliptic = [
                compensated_dot3([cosine, -sine, 0.0], [x, y, 0.0]),
                compensated_dot3([sine, cosine, 0.0], [x, y, 0.0]),
                z,
            ];
        }
        ecliptic
    } else if flags & FLG_NONUT != 0 {
        let m = if icrs { frames.p_nobias } else { frames.p };
        let geo_eq = mat3_vec3(m, geo)?;
        // Reduce longitude in the selected ecliptic coordinate system.
        let eps_mean = frames.eps_mean_rad.to_degrees().to_radians();
        rotate_equatorial_to_ecliptic(geo_eq, eps_mean)?
    } else if icrs {
        let geo_eq = mat3_vec3(frames.pn_nobias, geo)?;
        rotate_equatorial_to_ecliptic(geo_eq, frames.eps_true_rad)?
    } else {
        let geo_eq = mat3_vec3(frames.pn, geo)?;
        rotate_equatorial_to_ecliptic(geo_eq, frames.eps_true_rad)?
    };
    if want_xyz {
        Ok(v)
    } else {
        cartesian_to_spherical_deg(v)
    }
}

/// Convert WGS84 geodetic coordinates to ITRF Cartesian metres.
///
/// Longitude and latitude are radians; height is metres. Longitude is periodic
/// and is not clamped. Non-finite inputs are [`Error::NonFiniteArgument`].
pub fn geodetic_to_itrf_m(lon_rad: f64, lat_rad: f64, height_m: f64) -> Result<[f64; 3]> {
    geodetic_to_itrf_m_in(lon_rad, lat_rad, height_m, SiteEllipsoid::Wgs84)
}

/// Convert geodetic coordinates to ITRF metres on the selected ellipsoid.
///
/// The prime-vertical radius follows the ellipsoid's semimajor axis and
/// flattening. Longitude/latitude are radians and height is metres. WGS84
/// and IERS 2010 use their respective published ellipsoid parameters.
/// Non-finite inputs are [`Error::NonFiniteArgument`].
pub fn geodetic_to_itrf_m_in(
    lon_rad: f64,
    lat_rad: f64,
    height_m: f64,
    ellipsoid: SiteEllipsoid,
) -> Result<[f64; 3]> {
    for v in [lon_rad, lat_rad, height_m] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    Ok(geodetic_itrf_words(lon_rad, lat_rad, height_m, ellipsoid))
}

/// WGS84 geodetic radians/metres consumed by an actual array primitive.
/// Returns ITRF metres with the same ellipsoid arithmetic as the checked
/// geodetic API. Derived nonfinite operands propagate IEEE values; the array
/// owner supplies dtype, broadcasting and floating-point warning policy.
pub fn source_gd2gc_wgs84(lon_rad: f64, lat_rad: f64, height_m: f64) -> [f64; 3] {
    geodetic_itrf_words(lon_rad, lat_rad, height_m, SiteEllipsoid::Wgs84)
}

/// Shared ellipsoid formula after the consuming boundary's own admission.
/// Radian longitude/latitude and metre height produce ITRF metre components.
fn geodetic_itrf_words(
    lon_rad: f64,
    lat_rad: f64,
    height_m: f64,
    ellipsoid: SiteEllipsoid,
) -> [f64; 3] {
    let f = 1.0 / ellipsoid.inv_flattening();
    let e2 = 2.0 * f - f * f;
    let (slat, clat) = lat_rad.sin_cos();
    let (slon, clon) = lon_rad.sin_cos();
    let n = ellipsoid.semi_major_m() / (1.0 - e2 * slat * slat).sqrt();
    let r = (n + height_m) * clat;
    [r * clon, r * slon, ((1.0 - e2) * n + height_m) * slat]
}

/// IAU 2000 Earth Rotation Angle in radians, normalized to `[0, 2π)`.
///
/// The input is a finite UT1 Julian Day. The J2000 day offset and its
/// fractional rotation phase are kept separate to limit cancellation.
/// The defining IAU 2000 decimal constants are evaluated with compensated
/// products, retaining their phase at long epochs. The same rotation law is
/// implemented by ERFA `eraEra00`; attribution is in `docs/sources.md`.
/// Non-finite input is [`Error::NonFiniteArgument`].
pub fn earth_rotation_angle_rad(jd_ut1: f64) -> Result<f64> {
    if !jd_ut1.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_ut1 });
    }
    Ok(earth_rotation_turns_from_ut1_parts(jd_ut1, 0.0) * std::f64::consts::TAU)
}

/// Earth-rotation angle in radians at an acquired UT1 offset in days from
/// J2000. The split-date formula is shared with the checked scalar-JD entry;
/// no scalar JD is reconstructed and IEEE nonfinite operands are not refused.
pub fn source_earth_rotation_angle_from_j2000_offset(date2: f64) -> f64 {
    let days = date2.floor();
    earth_rotation_turns_from_split_days(days, date2 - days) * std::f64::consts::TAU
}

/// ERA in turns for two retained UT1 Julian-date parts, reduced to `[0,1)`.
///
/// Both operands are days. Separate their integral and fractional parts
/// before subtracting J2000, so sub-day corrections do not pass through a
/// large scalar epoch. Combining integral parts first also handles large
/// cancelling operands whose sum is an ordinary astronomical epoch. Integral
/// additions are exact within the supported astronomical domain. This helper
/// does not extend that physical domain; nonfinite/overflow arithmetic retains
/// IEEE propagation, while checked entry points validate their own operands.
pub(crate) fn earth_rotation_turns_from_ut1_parts(first: f64, second: f64) -> f64 {
    let first_day = first.floor();
    let second_day = second.floor();
    let days = (first_day + second_day) - FRAME_J2000;
    let fraction = (first - first_day) + (second - second_day);
    earth_rotation_turns_from_split_days(days, fraction)
}

/// Shared rotation law for an integral J2000 offset and its UT1 fraction.
/// Fractional parts may sum to a value in `[0,2]`, up to endpoint rounding;
/// unchecked nonfinite operands retain IEEE propagation. Return turns.
fn earth_rotation_turns_from_split_days(days: f64, fraction: f64) -> f64 {
    // These are the exact defining decimal constants minus their nearest
    // binary64 values, rounded once to binary64. They are not fitted phases.
    // The spin residual alone accumulates to 11.33 microarcseconds at 200 kyr.
    const SPIN_RESIDUAL: f64 = 1.1969606404306886e-19;
    const EPOCH_RESIDUAL: f64 = 1.6884986507648136e-17;

    let product = ERA_SPIN_TURNS_PER_DAY * days;
    // FMA recovers the product rounding before the large integer-turn part
    // is discarded. Reducing a rounded accumulated phase would lose up to
    // tens of microarcseconds even with a correctly specified UT1 clock.
    let product_error = ERA_SPIN_TURNS_PER_DAY.mul_add(days, -product);
    let correction = (SPIN_RESIDUAL.mul_add(days, product_error) + SPIN_RESIDUAL * fraction) % 1.0
        + EPOCH_RESIDUAL;
    let phase = (product % 1.0 + ERA_J2000_TURNS) + fraction + correction;
    // Keep the nominal one turn/day separate instead of rounding 1+spin.
    let phase = ERA_SPIN_TURNS_PER_DAY.mul_add(fraction, phase);
    let turns = phase.rem_euclid(1.0);
    // A negative remainder within half an ulp of zero can round its +1
    // adjustment to one. Canonicalize that endpoint on the same circle;
    // unchecked NaN propagation and every interior phase are unchanged.
    if turns == 1.0 { 0.0 } else { turns }
}

/// Active right-handed rotation about the z-axis by `angle_rad`.
pub(super) fn rot_z_rad(angle_rad: f64) -> [[f64; 3]; 3] {
    let (s, c) = angle_rad.sin_cos();
    [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]]
}

/// Polar-motion matrix with zero polar coordinates `xp = yp = 0`.
///
/// The IERS 2010 expression `W = R3(−s′)·R2(xp)·R1(yp)` then reduces to
/// `R3(−s′)`. `sp_rad` is the TIO locator in radians; a non-finite argument
/// is [`Error::NonFiniteArgument`].
pub fn polar_motion_matrix(sp_rad: f64) -> Result<[[f64; 3]; 3]> {
    if !sp_rad.is_finite() {
        return Err(Error::NonFiniteArgument { value: sp_rad });
    }
    Ok(rot_z_rad(-sp_rad))
}

/// Compose the IAU celestial-to-terrestrial rotation.
///
/// The polar-motion matrix multiplies Earth rotation and celestial pole
/// orientation: `c2t = rpom @ R3(era) @ rc2i`. Inputs are dimensionless matrices
/// and radians. All entries must be finite; matrix errors propagate. The
/// caller supplies a mutually consistent epoch, nutation and Earth rotation.
pub fn assemble_celestial_to_terrestrial(
    rpom: [[f64; 3]; 3],
    era_rad: f64,
    rc2i: [[f64; 3]; 3],
) -> Result<[[f64; 3]; 3]> {
    if !era_rad.is_finite() {
        return Err(Error::NonFiniteArgument { value: era_rad });
    }
    let theta = rot_z_rad(era_rad);
    mat3_mul(mat3_mul(rpom, theta)?, rc2i)
}

impl TopocentricOffset {
    /// Checked constructor: every component must be finite.
    pub fn new(pos_au: [f64; 3], vel_au_day: [f64; 3]) -> Result<Self> {
        check_finite_vec3(pos_au)?;
        check_finite_vec3(vel_au_day)?;
        Ok(Self { pos_au, vel_au_day })
    }

    /// Observer ICRS position offset from the geocenter, AU.
    pub fn pos_au(self) -> [f64; 3] {
        self.pos_au
    }

    /// Observer ICRS velocity offset from the geocenter velocity, AU/day.
    pub fn vel_au_day(self) -> [f64; 3] {
        self.vel_au_day
    }

    /// The `(position, velocity)` pair for the center-composition `topo` slot.
    pub fn pair(self) -> ([f64; 3], [f64; 3]) {
        (self.pos_au, self.vel_au_day)
    }
}

/// Geocentric observer position in ICRS AU and velocity in AU/day.
///
/// Geodetic longitude/latitude are degrees and height is metres. `c2t` is the
/// celestial-to-terrestrial matrix. The position is the transpose-rotated
/// ellipsoidal ITRF location; the velocity rotates the terrestrial `ω × r`
/// term with the same transpose. This nominal-rotation model omits polar-motion
/// rates and length-of-day corrections. WGS84 and IERS 2010 parameters are
/// selected by `ellipsoid`. Non-finite inputs or matrix values are errors.
pub fn topocentric_offset(
    lon_deg: f64,
    lat_deg: f64,
    height_m: f64,
    c2t: [[f64; 3]; 3],
    ellipsoid: SiteEllipsoid,
) -> Result<TopocentricOffset> {
    for v in [lon_deg, lat_deg, height_m] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    for row in c2t {
        check_finite_vec3(row)?;
    }
    let itrf = geodetic_to_itrf_m_in(
        lon_deg.to_radians(),
        lat_deg.to_radians(),
        height_m,
        ellipsoid,
    )?;
    topocentric_offset_from_itrf(itrf, c2t)
}

/// Rotate an already acquired ITRF position in metres into the ICRS site
/// position in AU and Earth-rotation velocity in AU/day. `c2t` is the owned
/// dimensionless celestial-to-terrestrial rotation. The left-associated dot,
/// AU conversion and angular-velocity tail are shared with
/// [`topocentric_offset`], allowing callers to acquire geodetic data before
/// their actual clock/frame operation. Nonfinite components are refused.
pub fn topocentric_offset_from_itrf(
    itrf: [f64; 3],
    c2t: [[f64; 3]; 3],
) -> Result<TopocentricOffset> {
    check_finite_vec3(itrf)?;
    for row in c2t {
        check_finite_vec3(row)?;
    }
    let (position, velocity) = source_topocentric_offset_from_itrf(itrf, c2t);
    TopocentricOffset::new(position, velocity)
}

/// Rotate an already acquired ITRF site into its geocentric celestial state.
///
/// Input position is metres and `c2t` is dimensionless. The output is ICRS AU
/// and nominal Earth-spin velocity in AU/day. This unchecked numerical helper
/// shares [`topocentric_offset_from_itrf`]'s geometry; non-finite components
/// propagate arithmetically rather than being validated here.
pub fn source_topocentric_offset_from_itrf(
    itrf: [f64; 3],
    c2t: [[f64; 3]; 3],
) -> ([f64; 3], [f64; 3]) {
    let (x_m, y_m, z_m) = (itrf[0], itrf[1], itrf[2]);
    // Transpose the celestial-to-terrestrial rotation and convert metres to AU.
    let t2c_0 = [c2t[0][0], c2t[1][0], c2t[2][0]];
    let t2c_1 = [c2t[0][1], c2t[1][1], c2t[2][1]];
    let t2c_2 = [c2t[0][2], c2t[1][2], c2t[2][2]];
    let pos_au = [
        (t2c_0[0] * x_m + t2c_0[1] * y_m + t2c_0[2] * z_m) / AU_IN_M,
        (t2c_1[0] * x_m + t2c_1[1] * y_m + t2c_1[2] * z_m) / AU_IN_M,
        (t2c_2[0] * x_m + t2c_2[1] * y_m + t2c_2[2] * z_m) / AU_IN_M,
    ];
    let vx_itrf = -EARTH_OMEGA_RAD_S * y_m;
    let vy_itrf = EARTH_OMEGA_RAD_S * x_m;
    let vz_itrf = 0.0;
    let vel_au_day = [
        (t2c_0[0] * vx_itrf + t2c_0[1] * vy_itrf + t2c_0[2] * vz_itrf) * SECONDS_PER_DAY / AU_IN_M,
        (t2c_1[0] * vx_itrf + t2c_1[1] * vy_itrf + t2c_1[2] * vz_itrf) * SECONDS_PER_DAY / AU_IN_M,
        (t2c_2[0] * vx_itrf + t2c_2[1] * vy_itrf + t2c_2[2] * vz_itrf) * SECONDS_PER_DAY / AU_IN_M,
    ];
    (pos_au, vel_au_day)
}

/// Whether a lunar point is defined geocentrically and has no observer
/// parallax: nodes, mean/osculating apogee and interpolated apsides.
pub fn is_topo_noop_point(body_id: i32) -> bool {
    matches!(
        body_id,
        MEAN_NODE | TRUE_NODE | MEAN_APOG | OSCU_APOG | INTP_APOG | INTP_PERG
    )
}

impl TopoPlan {
    /// The geodetic triple for [`TopoPlan::Topocentric`], else `None`.
    pub fn geopos(self) -> Option<(f64, f64, f64)> {
        match self {
            Self::Topocentric {
                lon_deg,
                lat_deg,
                height_m,
            } => Some((lon_deg, lat_deg, height_m)),
            Self::Geocentric | Self::TopoNoop => None,
        }
    }

    /// True only for [`TopoPlan::Topocentric`].
    pub fn uses_offset(self) -> bool {
        matches!(self, Self::Topocentric { .. })
    }
}

/// Resolve the topocentric flag and observer selection.
///
/// An explicit coordinate triple takes precedence over the context observer.
/// A real parallax request without either is [`Error::ConfigurationMissing`];
/// no geocentric substitution is made. Geocentrically defined lunar points
/// require no observer.
pub fn resolve_topo_request(
    flags: i32,
    body_id: i32,
    explicit_geopos: Option<(f64, f64, f64)>,
    observer: Option<Observer>,
) -> Result<TopoPlan> {
    if flags & FLG_TOPOCTR == 0 {
        return Ok(TopoPlan::Geocentric);
    }
    if is_topo_noop_point(body_id) {
        return Ok(TopoPlan::TopoNoop);
    }
    if let Some((lon_deg, lat_deg, height_m)) = explicit_geopos {
        return Ok(TopoPlan::Topocentric {
            lon_deg,
            lat_deg,
            height_m,
        });
    }
    if let Some(obs) = observer {
        return Ok(TopoPlan::Topocentric {
            lon_deg: obs.lon_deg(),
            lat_deg: obs.lat_deg(),
            height_m: obs.height_m(),
        });
    }
    Err(Error::ConfigurationMissing {
        missing_config: "observer_location",
        suggestion: "set an observer (set_topo) before requesting FLG_TOPOCTR",
    })
}

/// The frozen central `jd_tt − jd_ut1` difference (days) carried by the
/// speed-stencil layer.
///
/// Both entry shapes freeze `_dt` at the central epoch; every
/// `±_VEL_H` sample reuses it instead of re-converting scales.
pub fn frozen_delta_t_days(jd_tt: f64, jd_ut1: f64) -> Result<f64> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    if !jd_ut1.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_ut1 });
    }
    Ok(jd_tt - jd_ut1)
}

impl TopoSampler {
    /// Checked constructor: every field must be finite.
    pub fn new(
        lon_deg: f64,
        lat_deg: f64,
        height_m: f64,
        delta_t_days: f64,
        ellipsoid: SiteEllipsoid,
    ) -> Result<Self> {
        for v in [lon_deg, lat_deg, height_m, delta_t_days] {
            if !v.is_finite() {
                return Err(Error::NonFiniteArgument { value: v });
            }
        }
        Ok(Self {
            lon_deg,
            lat_deg,
            height_m,
            delta_t_days,
            ellipsoid,
        })
    }

    /// East longitude in degrees.
    pub fn lon_deg(self) -> f64 {
        self.lon_deg
    }

    /// Latitude in degrees.
    pub fn lat_deg(self) -> f64 {
        self.lat_deg
    }

    /// Height above the reference surface in metres.
    pub fn height_m(self) -> f64 {
        self.height_m
    }

    /// Frozen central `jd_tt − jd_ut1` in days.
    pub fn delta_t_days(self) -> f64 {
        self.delta_t_days
    }

    /// Reference ellipsoid for the site-to-ITRF conversion.
    pub fn ellipsoid(self) -> SiteEllipsoid {
        self.ellipsoid
    }

    /// UT1 epoch paired with a TT sample under the frozen difference.
    pub fn ut1_at(self, jd_tt_sample: f64) -> Result<f64> {
        if !jd_tt_sample.is_finite() {
            return Err(Error::NonFiniteArgument {
                value: jd_tt_sample,
            });
        }
        Ok(jd_tt_sample - self.delta_t_days)
    }

    /// Re-evaluate the offset at a TT sample epoch with the injected
    /// sample-epoch `c2t` (whose UT1 argument the entry layer derives via
    /// [`TopoSampler::ut1_at`], mirroring `topo_offset_fn`).
    pub fn offset_at(self, jd_tt_sample: f64, c2t: [[f64; 3]; 3]) -> Result<TopocentricOffset> {
        // Fail fast on a non-finite sample even though the UT1 value
        // travels with the matrix provider, not with this return.
        let _ = self.ut1_at(jd_tt_sample)?;
        topocentric_offset(
            self.lon_deg,
            self.lat_deg,
            self.height_m,
            c2t,
            self.ellipsoid,
        )
    }
}

/// Central-difference nutation-in-longitude rate, in degrees per TT day.
///
/// `dpsi_at` returns radians at a TT Julian Day. The fixed half-step is
/// [`VEL_H_DAYS`]; callers needing record-aware derivatives must select the
/// appropriate owned-domain sampler. Non-finite central epochs are errors;
/// callback failures propagate unchanged.
pub fn nutation_rate_deg_day(jd_tt: f64, dpsi_at: impl Fn(f64) -> Result<f64>) -> Result<f64> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let dpsi_m = dpsi_at(jd_tt - VEL_H_DAYS)?;
    let dpsi_p = dpsi_at(jd_tt + VEL_H_DAYS)?;
    Ok((dpsi_p - dpsi_m).to_degrees() / (2.0 * VEL_H_DAYS))
}

/// Linearly advance an observer-relative vector by `sample_offset_days`.
///
/// Position is AU and velocity is AU per TT day. This is a local tangent
/// approximation; downstream correction and frame boundaries validate their
/// inputs.
pub fn extrapolate_geometric_vector(
    geo_geom: [f64; 3],
    geo_vel: [f64; 3],
    sample_offset_days: f64,
) -> [f64; 3] {
    [
        geo_geom[0] + geo_vel[0] * sample_offset_days,
        geo_geom[1] + geo_vel[1] * sample_offset_days,
        geo_geom[2] + geo_vel[2] * sample_offset_days,
    ]
}

/// Replace a linear site extrapolation by the actual sampled site.
///
/// Subtracts `topo_sample − topo_central − topo_central_vel·s`, in AU, from
/// the extrapolated relative vector. The caller supplies both site positions
/// and the central AU/day velocity for the same offset `s` in days.
pub fn remove_topo_curvature(
    extrapolated: [f64; 3],
    topo_central_pos: [f64; 3],
    topo_central_vel: [f64; 3],
    topo_sample_pos: [f64; 3],
    sample_offset_days: f64,
) -> [f64; 3] {
    [
        extrapolated[0]
            - (topo_sample_pos[0] - topo_central_pos[0] - topo_central_vel[0] * sample_offset_days),
        extrapolated[1]
            - (topo_sample_pos[1] - topo_central_pos[1] - topo_central_vel[1] * sample_offset_days),
        extrapolated[2]
            - (topo_sample_pos[2] - topo_central_pos[2] - topo_central_vel[2] * sample_offset_days),
    ]
}

/// One-way light time in days from the sampled vector's Euclidean AU norm.
/// Recomputing it keeps the aberration normalization consistent with the ray.
pub fn sample_light_time_days(extrapolated: [f64; 3]) -> f64 {
    vec3_dist(extrapolated) / C_LIGHT_AU_DAY
}

impl StencilTopo {
    /// Checked constructor: every supplied component must be finite.
    pub fn new(
        central_pos: [f64; 3],
        central_vel: [f64; 3],
        sample: Option<([f64; 3], [f64; 3])>,
    ) -> Result<Self> {
        check_finite_vec3(central_pos)?;
        check_finite_vec3(central_vel)?;
        if let Some((pos, vel)) = sample {
            check_finite_vec3(pos)?;
            check_finite_vec3(vel)?;
        }
        Ok(Self {
            central_pos,
            central_vel,
            sample,
        })
    }

    /// Central-epoch topocentric position offset, ICRS AU.
    pub fn central_pos(self) -> [f64; 3] {
        self.central_pos
    }

    /// Central-epoch topocentric velocity offset, ICRS AU/day.
    pub fn central_vel(self) -> [f64; 3] {
        self.central_vel
    }

    /// Sample-epoch `(position, velocity)` pair, if re-evaluated.
    pub fn sample(self) -> Option<([f64; 3], [f64; 3])> {
        self.sample
    }
}

impl StencilCenter {
    /// Checked constructor: every vector component must be finite.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        geo_geom: [f64; 3],
        geo_vel: [f64; 3],
        observer: [f64; 3],
        observer_vel: [f64; 3],
        defl_delta: [f64; 3],
        target_body_id: i32,
        do_defl: bool,
        do_aber: bool,
        topo: Option<StencilTopo>,
    ) -> Result<Self> {
        check_finite_vec3(geo_geom)?;
        check_finite_vec3(geo_vel)?;
        check_finite_vec3(observer)?;
        check_finite_vec3(observer_vel)?;
        check_finite_vec3(defl_delta)?;
        Ok(Self {
            geo_geom,
            geo_vel,
            observer,
            observer_vel,
            defl_delta,
            target_body_id,
            do_defl,
            do_aber,
            topo,
        })
    }

    /// Pre-deflection, pre-aberration apparent vector at `jd_tt`, AU.
    pub fn geo_geom(self) -> [f64; 3] {
        self.geo_geom
    }

    /// Analytic light-time-corrected rate, AU/day.
    pub fn geo_vel(self) -> [f64; 3] {
        self.geo_vel
    }

    /// Already-selected observer position, AU.
    pub fn observer(self) -> [f64; 3] {
        self.observer
    }

    /// Already-selected observer velocity, AU/day.
    pub fn observer_vel(self) -> [f64; 3] {
        self.observer_vel
    }

    /// Captured position-path deflection increment, AU.
    pub fn defl_delta(self) -> [f64; 3] {
        self.defl_delta
    }

    /// Target body id selecting the Sun freeze branch.
    pub fn target_body_id(self) -> i32 {
        self.target_body_id
    }

    /// Whether the stencil re-applies gravitational deflection.
    pub fn do_defl(self) -> bool {
        self.do_defl
    }

    /// Whether the stencil re-applies aberration.
    pub fn do_aber(self) -> bool {
        self.do_aber
    }

    /// Central-plus-sample topocentric state, if topocentric.
    pub fn topo(self) -> Option<StencilTopo> {
        self.topo
    }
}

#[cfg(test)]
mod source_observer_tests {
    use super::*;

    #[test]
    fn acquired_observer_words_share_checked_formula_without_its_admission() {
        let xyz = source_gd2gc_wgs84(0.25, -0.75, 12.);
        assert_eq!(xyz, geodetic_to_itrf_m(0.25, -0.75, 12.).unwrap());
        let matrix = crate::erfa::c2t06a(2451545., 2451545., 0., 0.).unwrap();
        let raw = source_topocentric_offset_from_itrf(xyz, matrix);
        assert_eq!(
            raw,
            topocentric_offset_from_itrf(xyz, matrix).unwrap().pair()
        );
        assert!(source_gd2gc_wgs84(f64::NAN, 0., 0.)[0].is_nan());
        assert!(geodetic_to_itrf_m(f64::NAN, 0., 0.).is_err());
        assert!(source_topocentric_offset_from_itrf([f64::NAN; 3], matrix).0[0].is_nan());
        assert!(topocentric_offset_from_itrf([f64::NAN; 3], matrix).is_err());
        assert_eq!(
            source_earth_rotation_angle_from_j2000_offset(0.),
            earth_rotation_angle_rad(FRAME_J2000).unwrap()
        );
    }
}

#[cfg(test)]
mod compensated_product_tests {
    use super::*;

    #[test]
    fn shared_frame_components_preserve_every_binary64_word() {
        for jd in [
            -3_100_000.5,
            -105149.5,
            1721060.5,
            2451545.0,
            2634165.5,
            8000014.5,
        ] {
            for (psi, deps) in [(0.0, 0.0), (-8e-5, 4e-5), (1e-4, -2e-5)] {
                let frames = frame_matrices_for_date(jd, psi, deps).unwrap();
                let bits = |m: [[f64; 3]; 3]| m.map(|row| row.map(f64::to_bits));
                for (bias, precession, combined) in [
                    (true, frames.p(), frames.pn()),
                    (false, frames.p_nobias(), frames.pn_nobias()),
                ] {
                    let (expected, obliquity) = vondrak_pn_matrix(jd, psi, deps, bias).unwrap();
                    assert_eq!(bits(combined), bits(expected));
                    assert_eq!(
                        bits(precession),
                        bits(vondrak_precession_matrix(jd, bias).unwrap())
                    );
                    assert_eq!(frames.eps_true_rad().to_bits(), obliquity.to_bits());
                }
            }
        }
        for (jd, psi, deps, invalid) in [
            (f64::INFINITY, f64::NEG_INFINITY, 0.0, f64::INFINITY),
            (
                2451545.0,
                f64::NEG_INFINITY,
                f64::INFINITY,
                f64::NEG_INFINITY,
            ),
            (2451545.0, 0.0, f64::INFINITY, f64::INFINITY),
        ] {
            assert!(matches!(frame_matrices_for_date(jd, psi, deps),
                Err(Error::NonFiniteArgument { value }) if value == invalid));
        }
    }
    #[test]
    fn matrix_and_vector_products_retain_an_exact_cancellation_residual() {
        let epsilon = 2.0_f64.powi(-27);
        let row = [-1.0, 1.0 + epsilon, 0.0];
        let column = [1.0, 1.0 - epsilon, 0.0];
        let exact = -2.0_f64.powi(-54);
        let matrix = [row, [0.0; 3], [0.0; 3]];
        assert_eq!(mat3_vec3(matrix, column).unwrap()[0], exact);
        assert_eq!(
            mat3_mul(
                matrix,
                [[1.0, 0.0, 0.0], [1.0 - epsilon, 0.0, 0.0], [0.0; 3]]
            )
            .unwrap()[0][0],
            exact
        );
    }

    #[test]
    fn native_ecliptic_latitude_is_independent_of_nutation() {
        let ray = [
            0.17714781512074745,
            -0.8874263370653095,
            -0.3847419511610094,
        ];
        for epoch in [-105149.5, 1721060.5, 2451545.000740741, 2634165.5] {
            for icrs in [0, FLG_ICRS] {
                let frames = frame_matrices_for_date(epoch, -8e-5, 4e-5).unwrap();
                let mean = frame_transform(ray, epoch, icrs | FLG_NONUT, true, frames).unwrap();
                let true_place = frame_transform(ray, epoch, icrs, true, frames).unwrap();
                assert_eq!(
                    true_place[2], mean[2],
                    "nutation cannot change the ecliptic pole"
                );
                let reference = if icrs == 0 {
                    mat3_vec3(FRAME_BIAS_LTPB, ray).unwrap()
                } else {
                    ray
                };
                assert!(
                    (vec3_dist(true_place) - vec3_dist(reference)).abs()
                        <= 32.0 * f64::EPSILON * vec3_dist(reference)
                );
                let expected_z = compensated_dot3(ecliptic_pole_vondrak(epoch), reference);
                assert_eq!(true_place[2], expected_z);
                // Equatorial consumers retain the actual supplied PN matrix.
                assert_eq!(
                    frame_transform(ray, epoch, FLG_EQUATORIAL, true, frames).unwrap(),
                    mat3_vec3(frames.pn(), ray).unwrap()
                );
            }
        }
    }

    #[test]
    fn custom_frame_matrices_retain_their_explicit_coordinate_system() {
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let frames =
            FrameMatrices::new(identity, identity, identity, identity, 0.41, 0.409).unwrap();
        let ray = [0.25, 0.5, 1.0];
        assert_eq!(
            frame_transform(ray, FRAME_J2000, 0, true, frames).unwrap(),
            rotate_equatorial_to_ecliptic(ray, 0.41).unwrap()
        );
    }
}

#[cfg(test)]
mod short_rotation_interval_tests {
    use super::*;

    #[test]
    fn earth_rotation_resolves_a_short_interval_after_many_whole_turns() {
        // The IAU law has the exact constant derivative 2π(1 + spin).
        // A dyadic displacement is representable at each tested epoch.
        let step = 2.0_f64.powi(-20);
        for day in [-1_000_000.0, -18263.0, 0.0, 18263.0, 1_000_000.0] {
            let epoch = day + 0.25;
            let previous = source_earth_rotation_angle_from_j2000_offset(epoch - step);
            let following = source_earth_rotation_angle_from_j2000_offset(epoch + step);
            let expected = 2.0 * step * std::f64::consts::TAU * (1.0 + ERA_SPIN_TURNS_PER_DAY);
            let delta = (following - previous + std::f64::consts::PI)
                .rem_euclid(std::f64::consts::TAU)
                - std::f64::consts::PI;
            assert!(
                (delta - expected).abs() <= 8.0 * f64::EPSILON,
                "day {day}: rotation {delta}, exact {expected}"
            );
        }
    }

    #[test]
    fn reduced_rotation_keeps_the_j2000_phase_and_periodic_range() {
        assert_eq!(
            source_earth_rotation_angle_from_j2000_offset(0.0),
            ERA_J2000_TURNS * std::f64::consts::TAU
        );
        for epoch in [-1_000_000.75, -1.0, -0.25, 0.25, 1.0, 1_000_000.75] {
            let angle = source_earth_rotation_angle_from_j2000_offset(epoch);
            assert!((0.0..std::f64::consts::TAU).contains(&angle));
        }
    }
}

#[cfg(test)]
mod exact_frame_cache_tests {
    use super::*;

    fn words(frames: FrameMatrices) -> Vec<u64> {
        let mut words: Vec<_> = [frames.p, frames.p_nobias, frames.pn, frames.pn_nobias]
            .into_iter()
            .flatten()
            .flatten()
            .map(f64::to_bits)
            .collect();
        words.extend([frames.eps_mean_rad.to_bits(), frames.eps_true_rad.to_bits()]);
        if let Some((matrix, psi)) = frames.native_ecliptic {
            words.extend(matrix.into_iter().flatten().map(f64::to_bits));
            words.push(psi.to_bits());
        }
        words
    }

    #[test]
    fn every_operand_and_thread_keeps_the_uncached_frame_words() {
        std::thread::scope(|scope| {
            for thread in 0..8 {
                scope.spawn(move || {
                    for epoch in [1000000.25, 2415020.5, 2451545.0, 5000000.75] {
                        for psi in [0.0, -0.0, 0.00001, 0.000010000000000000002] {
                            let deps = thread as f64 * 1e-6;
                            let expected =
                                words(frame_matrices_uncached(epoch, psi, deps).unwrap());
                            for _ in 0..3 {
                                assert_eq!(
                                    words(frame_matrices_for_date(epoch, psi, deps).unwrap()),
                                    expected
                                );
                            }
                        }
                    }
                });
            }
        });
        for (tt, psi, deps) in [
            (f64::NAN, 0.0, 0.0),
            (2451545.0, f64::INFINITY, 0.0),
            (2451545.0, 0.0, f64::NAN),
        ] {
            assert!(frame_matrices_for_date(tt, psi, deps).is_err());
        }
    }
}

#[cfg(test)]
mod epoch_pair_rate_tests {
    use super::*;

    #[test]
    fn complete_epoch_pair_rates_match_independent_mp120_derivatives() {
        // Full 8/14-term MP pole series and atan2 Jacobian, not the general
        // precession polynomial. Every distinct clock/law from the complete
        // 197-request scalar-offset qualification is represented here.
        let cases: &[(f64, f64, f64, bool, f64)] = &[
            (
                2451545.000740741,
                24.042044444444457,
                2433282.42345905,
                false,
                3.824470779726048e-05,
            ),
            (
                2451545.0,
                24.042044444444457,
                2433282.42345905,
                false,
                3.824470779691945e-05,
            ),
            (
                2451545.000740741,
                23.85709166666667,
                2451545.0,
                false,
                3.824470441615096e-05,
            ),
            (
                2451545.000740741,
                26.413055555555555,
                2415020.31352,
                false,
                3.824471069157877e-05,
            ),
            (
                2451545.000740741,
                21.014166666666668,
                2415020.0,
                false,
                3.8244710691624266e-05,
            ),
            (
                2451545.000740741,
                0.0,
                1925312.0,
                false,
                3.824460046728006e-05,
            ),
            (
                2451545.000740741,
                22.36388888888889,
                2415020.0,
                false,
                3.8244710691624266e-05,
            ),
            (
                2451545.000740741,
                30.0,
                2494279.0,
                false,
                3.8244694614950074e-05,
            ),
            (
                2451545.000740741,
                20.91,
                2412908.1187773715,
                false,
                3.824471099483657e-05,
            ),
            (
                2451545.000740741,
                0.0,
                1854089.0,
                false,
                3.824455345821144e-05,
            ),
            (
                2451545.000740741,
                0.0,
                1834000.0,
                false,
                3.824453875387688e-05,
            ),
            (
                2451545.000740741,
                0.0,
                1797110.0,
                false,
                3.824451009422242e-05,
            ),
            (
                2451545.000740741,
                0.0,
                1774464.0,
                false,
                3.824449143705729e-05,
            ),
            (
                2451545.000740741,
                -4.466666666666667,
                1684613.0,
                false,
                3.8244409439999925e-05,
            ),
            (
                2451545.000740741,
                0.0,
                1807702.0,
                false,
                3.824451854281564e-05,
            ),
            (
                2451545.000740741,
                0.0,
                1920199.0,
                false,
                3.824459735856549e-05,
            ),
            (
                2451545.000740741,
                0.0,
                1927136.0,
                false,
                3.824460156632965e-05,
            ),
            (
                2451545.000740741,
                0.0,
                1903395.1338202194,
                false,
                3.824458685216622e-05,
            ),
            (
                2451545.000740741,
                -0.21510829031467438,
                1903396.8775657457,
                false,
                3.824458685327952e-05,
            ),
            (
                2451545.000740741,
                -0.23810671619139612,
                1903396.8542324193,
                false,
                3.824458685326463e-05,
            ),
            (
                2451545.000740741,
                -0.7922573227712064,
                1903395.1338202194,
                false,
                3.824458685216622e-05,
            ),
            (
                2451545.000740741,
                2.1105983371992068,
                1903395.1338202194,
                false,
                3.824458685216622e-05,
            ),
            (
                2451545.000740741,
                0.0,
                1911798.0,
                false,
                3.82445921614738e-05,
            ),
            (
                2451545.000740741,
                0.0,
                1805145.0,
                false,
                3.824451651946058e-05,
            ),
            (
                2451545.000740741,
                25.0,
                2451544.5,
                false,
                3.824470441625018e-05,
            ),
            (
                2451545.000740741,
                -2.9422,
                1775845.5007407407,
                false,
                3.824449259839909e-05,
            ),
            (
                2451545.000740741,
                0.0,
                1825599.0,
                false,
                3.824453241595975e-05,
            ),
            (
                2451545.000740741,
                0.0,
                1825233.3060523525,
                false,
                3.824453213754304e-05,
            ),
            (
                2451545.000740741,
                0.0,
                1827424.7590899938,
                false,
                3.824453380282397e-05,
            ),
            (
                2451545.000740741,
                23.245339462821466,
                2435553.5,
                false,
                3.824470740325233e-05,
            ),
            (
                2451545.000740741,
                24.0,
                2451545.0,
                false,
                3.824470441615096e-05,
            ),
            (
                2460310.5,
                23.85709166666667,
                2451545.0,
                false,
                3.824874036670925e-05,
            ),
            (
                2460310.5,
                23.85709166666667,
                2451545.0,
                true,
                3.824873047202019e-05,
            ),
            (
                2451545.0,
                23.85709166666667,
                2451545.0,
                false,
                3.82447044158099e-05,
            ),
            (
                2460310.500740741,
                23.85709166666667,
                2451545.0,
                false,
                3.824874036705031e-05,
            ),
        ];
        for &(tt, value, epoch, ecl_date, expected) in cases {
            let actual = epoch_pair_ayanamsha_rate_deg_day(tt, value, epoch, ecl_date).unwrap();
            assert!(
                (actual - expected).abs() <= 32.0 * f64::EPSILON * expected.abs(),
                "tt={tt}, epoch={epoch}, ecl_date={ecl_date}: {actual} != {expected}"
            );
        }
    }

    #[test]
    fn epoch_pair_rates_reject_nonfinite_selected_operands() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            for (tt, value, epoch) in [
                (bad, 24.0, 2451545.0),
                (2451545.0, bad, 2451545.0),
                (2451545.0, 24.0, bad),
            ] {
                assert!(matches!(
                    epoch_pair_ayanamsha_rate_deg_day(tt, value, epoch, false),
                    Err(Error::NonFiniteArgument { .. })
                ));
            }
        }
    }
}
