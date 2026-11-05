//! Equinox-based precession, nutation and terrestrial observer frames.
//!
//! The numerical models are IAU 2000A nutation, Capitaine P03 precession,
//! IERS 2010 complementary terms, equinox-method sidereal time and WGS84
//! geodetic coordinates. Barycentric positions and velocities are in AU and
//! AU/day. Angles are radians unless a function explicitly names another unit.
//! Split TT, TDB and UT1 epochs retain the precision of short time intervals.
//! This chain excludes polar motion; callers requiring it must use the IERS
//! celestial-to-terrestrial model.

use crate::errors::Error;
use crate::iau2000_tables::{
    IAU2000_KE0_T, IAU2000_KE1, IAU2000_LUNISOLAR_LON, IAU2000_LUNISOLAR_OBL, IAU2000_NALS_T,
    IAU2000_NAPL_T, IAU2000_NUT_LON, IAU2000_NUT_OBL, IAU2000_SE0_T0, IAU2000_SE0_T1,
};

mod velocity;
pub(crate) use velocity::iau2000a_rate_tenth_microarcsec_per_day;

/// Arcseconds to radians.
pub const FRAME_ASEC2RAD: f64 = 4.84813681109536e-6;
/// Arcseconds per full circle.
pub const FRAME_ASEC360: f64 = 1296000.0;
/// J2000 Julian date.
pub const FRAME_EPOCH_J2000: f64 = 2451545.0;
/// Seconds per day.
pub const FRAME_DAY_S: f64 = 86400.0;
/// Metres per astronomical unit, IAU 2012.
pub const FRAME_AU_M: f64 = 149597870700.0;
/// Speed of light in metres per second, exact SI value.
pub const FRAME_C_M_S: f64 = 299792458.0;
/// Speed of light in astronomical units per day.
pub const FRAME_LIGHT_SPEED_AU_DAY_LOCAL: f64 = FRAME_C_M_S * FRAME_DAY_S / FRAME_AU_M;
/// Nominal Earth angular velocity in radians per second.
pub const FRAME_ANGVEL_RAD_S: f64 = 7.2921150e-5;
/// Mean Earth radius in metres.
pub const FRAME_EARTH_RADIUS_M: f64 = 6378136.6;
/// Solar gravitational parameter in m³/s², DE405.
pub const DE405_SOLAR_GM_M3_S2: f64 = 1.32712440017987e+20;
/// Tenth-microarcsecond to radians.
pub const FRAME_TENTH_USEC_2_RAD: f64 = FRAME_ASEC2RAD / 1e7;
/// Radians per full circle.
pub const FRAME_TAU: f64 = std::f64::consts::TAU;
/// WGS84 semimajor axis in metres.
pub const FRAME_WGS84_A_M: f64 = 6378137.0;
/// WGS84 inverse flattening.
pub const FRAME_WGS84_INV_F: f64 = 298.257223563;
/// Linear-in-century complementary sine-series coefficient in arcseconds.
pub const IAU2000_SE1_0: f64 = -0.87e-6;
/// Linear-in-century complementary cosine-series coefficient in arcseconds.
pub const IAU2000_SE1_1: f64 = 0.0;

/// Fundamental-argument polynomials for l, l′, F, D and Ω.
/// Each row contains constant through fourth-degree coefficients in arcseconds.
pub(super) const FUND_ROWS: [[f64; 5]; 5] = [
    [
        485868.249036,
        1717915923.2178,
        31.8792,
        0.051635,
        -0.00024470,
    ],
    [
        1287104.79305,
        129596581.0481,
        -0.5532,
        0.000136,
        -0.00001149,
    ],
    [
        335779.526232,
        1739527262.8478,
        -12.7512,
        -0.001037,
        0.00000417,
    ],
    [
        1072260.70369,
        1602961601.2090,
        -6.3706,
        0.006593,
        -0.00003169,
    ],
    [450160.398036, -6962890.5431, 7.4722, 0.007702, -0.00005939],
];

/// Planetary-argument constants and rates in radians and radians per century.
pub(super) const ANOMALY_ROWS: [[f64; 2]; 14] = [
    [2.35555598, 8328.6914269554],
    [6.24006013, 628.301955],
    [1.627905234, 8433.466158131],
    [5.198466741, 7771.3771468121],
    [2.18243920, -33.757045],
    [4.402608842, 2608.7903141574],
    [3.176146697, 1021.3285546211],
    [1.753470314, 628.3075849991],
    [6.203480913, 334.0612426700],
    [0.599546497, 52.9690962641],
    [0.874016757, 21.3299104960],
    [5.481293871, 7.4781598567],
    [5.321159000, 3.8127774000],
    [0.02438175, 0.00000538691],
];

/// Sine of an angle in radians. Non-finite inputs yield NaN.
pub(crate) fn scalar_sin(angle: f64) -> f64 {
    angle.sin()
}

/// Cosine of an angle in radians. Non-finite inputs yield NaN.
pub(crate) fn scalar_cos(angle: f64) -> f64 {
    angle.cos()
}

/// Euclidean remainder for a positive divisor, used for angular reductions.
pub fn floored_mod(a: f64, b: f64) -> f64 {
    a.rem_euclid(b)
}

/// Signed remainder for a positive divisor, used by fundamental arguments.
pub fn truncated_mod(a: f64, b: f64) -> f64 {
    a % b
}

fn check_finite(value: f64) -> Result<f64, Error> {
    if !value.is_finite() {
        return Err(Error::NonFiniteArgument { value });
    }
    Ok(value)
}

/// USNO Circular 179 equation 2.6, TDB minus TT in seconds.
/// The epoch is a two-part Julian date. Non-finite inputs return an error.
/// The shared time-scale implementation supplies the scalar series.
pub fn tdb_minus_tt_seconds(jd_whole: f64, fraction: f64) -> Result<f64, Error> {
    crate::time::tdb_minus_tt_seconds(jd_whole, fraction)
}

/// Capitaine et al. (2003), equation 39, mean ecliptic obliquity in arcseconds.
/// The input is a finite TDB Julian date; non-finite dates return an error.
pub fn mean_obliquity_arcsec(jd_tdb: f64) -> Result<f64, Error> {
    let jd_tdb = check_finite(jd_tdb)?;
    let delta = jd_tdb - FRAME_EPOCH_J2000;
    let t = delta / 36525.0;
    let value = -0.0000000434 * t;
    let value = value - 0.000000576;
    let value = value * t;
    let value = value + 0.00200340;
    let value = value * t;
    let value = value - 0.0001831;
    let value = value * t;
    let value = value - 46.836769;
    let value = value * t;
    Ok(value + 84381.406)
}

/// Delaunay arguments l, l′, F, D and Ω in radians.
/// The input is finite Julian centuries from J2000; signed full-turn reduction
/// is applied after each polynomial. Non-finite epochs return an error.
pub fn fundamental_arguments_radians(t_centuries: f64) -> Result<[f64; 5], Error> {
    let t = check_finite(t_centuries)?;
    let mut out = [0.0; 5];
    for (row, slot) in FUND_ROWS.iter().zip(out.iter_mut()) {
        // Horner evaluation of each Delaunay-argument polynomial.
        let mut a = row[4] * t;
        a += row[3];
        a *= t;
        a += row[2];
        a *= t;
        a += row[1];
        a *= t;
        a += row[0];
        *slot = truncated_mod(a, FRAME_ASEC360) * FRAME_ASEC2RAD;
    }
    Ok(out)
}

/// Dot product of equal-length dimensionless series.
///
/// Finite inputs are accumulated in their original order; arithmetic overflow
/// is not checked, so a successful result need not be finite.
///
/// # Errors
///
/// Returns [`Error::NonFiniteArgument`] for a consumed nonfinite operand.
///
/// # Panics
///
/// Panics if `x` and `y` have different lengths.
pub fn series_dot(x: &[f64], y: &[f64]) -> Result<f64, Error> {
    assert_eq!(x.len(), y.len());
    let mut sum = 0.0;
    for (&a, &b) in x.iter().zip(y) {
        check_finite(a)?;
        check_finite(b)?;
        sum += a * b;
    }
    Ok(sum)
}

/// Dot product of a coefficient-table column and a vector.
///
/// Uses column `col` from each row and the corresponding value in `x`.
/// Only consumed column values are checked; arithmetic overflow is not checked.
///
/// # Errors
///
/// Returns [`Error::NonFiniteArgument`] for a consumed nonfinite operand.
///
/// # Panics
///
/// Panics if the table and vector lengths differ or `col >= C`.
pub fn strided_col_dot<const C: usize>(
    table: &[[f64; C]],
    col: usize,
    x: &[f64],
) -> Result<f64, Error> {
    assert_eq!(table.len(), x.len());
    assert!(col < C);
    let mut sum = 0.0;
    for (row, &value) in table.iter().zip(x) {
        check_finite(value)?;
        check_finite(row[col])?;
        sum += value * row[col];
    }
    Ok(sum)
}

/// Integer coefficient matrix times a dimensionless vector in binary64.
///
/// Each row must have `cols` entries and the iterator must yield exactly
/// `rows` rows. Finite weights do not guarantee finite products or sums.
///
/// # Errors
///
/// Returns [`Error::NonFiniteArgument`] for a nonfinite weight.
///
/// # Panics
///
/// Panics if `weights.len() != cols`, a row has the wrong length, or the
/// iterator yields a different number of rows than `rows`.
pub fn integer_matrix_vector(
    matrix: impl IntoIterator<Item = impl AsRef<[i32]>>,
    rows: usize,
    cols: usize,
    weights: &[f64],
) -> Result<Vec<f64>, Error> {
    assert_eq!(weights.len(), cols);
    for &value in weights {
        check_finite(value)?;
    }
    let mut out = Vec::with_capacity(rows);
    for row in matrix {
        let row = row.as_ref();
        assert_eq!(row.len(), cols);
        out.push(
            row.iter()
                .zip(weights)
                .map(|(&a, &b)| f64::from(a) * b)
                .sum(),
        );
    }
    assert_eq!(out.len(), rows);
    Ok(out)
}

/// Row-major binary64 matrix times a dimensionless vector.
///
/// `flat` contains `rows * cols` coefficients in row-major order. Each
/// consumed row is reduced by [`series_dot`]; overflow is not checked.
///
/// # Errors
///
/// Returns [`Error::NonFiniteArgument`] for a consumed nonfinite coefficient
/// or weight. With zero rows no operands are consumed.
///
/// # Panics
///
/// Panics if the flat matrix length differs from `rows * cols` or the weight
/// length differs from `cols`. The dimension product can itself panic on
/// integer overflow when overflow checks are enabled.
pub fn matrix_vector(
    flat: &[f64],
    rows: usize,
    cols: usize,
    weights: &[f64],
) -> Result<Vec<f64>, Error> {
    assert_eq!(flat.len(), rows * cols);
    assert_eq!(weights.len(), cols);
    (0..rows)
        .map(|row| series_dot(&flat[row * cols..(row + 1) * cols], weights))
        .collect()
}

/// IAU 2000A nutation in longitude and obliquity, in tenth-microarcseconds.
/// The input is a finite TT Julian date. Luni-solar and planetary terms use
/// the bundled IAU coefficients; non-finite dates return an error.
pub fn iau2000a_tenth_microarcsec(jd_tt: f64) -> Result<(f64, f64), Error> {
    let jd_tt = check_finite(jd_tt)?;
    if let Some(angles) = crate::pure_cache::reference_nutation(jd_tt.to_bits()) {
        return Ok(angles);
    }
    let angles = iau2000a_uncached(jd_tt)?;
    crate::pure_cache::retain_reference_nutation(jd_tt.to_bits(), angles);
    Ok(angles)
}

fn iau2000a_uncached(jd_tt: f64) -> Result<(f64, f64), Error> {
    let t = (jd_tt - FRAME_EPOCH_J2000) / 36525.0;
    let fund = fundamental_arguments_radians(t)?;
    // Luni-solar argument matrix (678 effective rows).
    let arg = integer_matrix_vector(IAU2000_NALS_T.iter(), IAU2000_NALS_T.len(), 5, &fund)?;
    let mut sarg = Vec::with_capacity(arg.len());
    let mut carg = Vec::with_capacity(arg.len());
    for &a in arg.iter() {
        sarg.push(scalar_sin(a));
        carg.push(scalar_cos(a));
    }
    let mut dpsi = strided_col_dot(&IAU2000_LUNISOLAR_LON, 0, &sarg)?
        + strided_col_dot(&IAU2000_LUNISOLAR_LON, 1, &sarg)? * t
        + strided_col_dot(&IAU2000_LUNISOLAR_LON, 2, &carg)?;
    let mut deps = strided_col_dot(&IAU2000_LUNISOLAR_OBL, 0, &carg)?
        + strided_col_dot(&IAU2000_LUNISOLAR_OBL, 1, &carg)? * t
        + strided_col_dot(&IAU2000_LUNISOLAR_OBL, 2, &sarg)?;
    // Planetary arguments: `t * rate + constant`, last row times `t`.
    let mut parg = [0.0; 14];
    for (i, row) in ANOMALY_ROWS.iter().enumerate() {
        parg[i] = t * row[1] + row[0];
    }
    parg[13] *= t;
    let parg_vec = integer_matrix_vector(IAU2000_NAPL_T.iter(), IAU2000_NAPL_T.len(), 14, &parg)?;
    let mut sparg = Vec::with_capacity(parg_vec.len());
    let mut carpg = Vec::with_capacity(parg_vec.len());
    for &a in parg_vec.iter() {
        sparg.push(scalar_sin(a));
        carpg.push(scalar_cos(a));
    }
    dpsi += strided_col_dot(&IAU2000_NUT_LON, 0, &sparg)?
        + strided_col_dot(&IAU2000_NUT_LON, 1, &carpg)?;
    deps += strided_col_dot(&IAU2000_NUT_OBL, 0, &sparg)?
        + strided_col_dot(&IAU2000_NUT_OBL, 1, &carpg)?;
    Ok((dpsi, deps))
}

/// IERS complementary equation-of-the-equinoxes terms in radians.
/// The input is a finite TT Julian date; non-finite dates return an error.
pub fn equinoxes_complimentary_terms_radians(jd_tt: f64) -> Result<f64, Error> {
    let jd_tt = check_finite(jd_tt)?;
    let t = (jd_tt - FRAME_EPOCH_J2000) / 36525.0;
    let mut fa = [0.0; 14];
    fa[0] = (485868.249036
        + (715923.2178 + (31.8792 + (0.051635 + (-0.00024470) * t) * t) * t) * t)
        * FRAME_ASEC2RAD
        + floored_mod(1325.0 * t, 1.0) * FRAME_TAU;
    fa[1] = (1287104.793048
        + (1292581.0481 + (-0.5532 + (0.000136 + (-0.00001149) * t) * t) * t) * t)
        * FRAME_ASEC2RAD
        + floored_mod(99.0 * t, 1.0) * FRAME_TAU;
    fa[2] = (335779.526232
        + (295262.8478 + (-12.7512 + (-0.001037 + (0.00000417) * t) * t) * t) * t)
        * FRAME_ASEC2RAD
        + floored_mod(1342.0 * t, 1.0) * FRAME_TAU;
    fa[3] = (1072260.703692
        + (1105601.2090 + (-6.3706 + (0.006593 + (-0.00003169) * t) * t) * t) * t)
        * FRAME_ASEC2RAD
        + floored_mod(1236.0 * t, 1.0) * FRAME_TAU;
    fa[4] = (450160.398036
        + (-482890.5431 + (7.4722 + (0.007702 + (-0.00005939) * t) * t) * t) * t)
        * FRAME_ASEC2RAD
        + floored_mod(-5.0 * t, 1.0) * FRAME_TAU;
    fa[5] = 4.402608842 + 2608.7903141574 * t;
    fa[6] = 3.176146697 + 1021.3285546211 * t;
    fa[7] = 1.753470314 + 628.3075849991 * t;
    fa[8] = 6.203480913 + 334.0612426700 * t;
    fa[9] = 0.599546497 + 52.9690962641 * t;
    fa[10] = 0.874016757 + 21.3299104960 * t;
    fa[11] = 5.481293872 + 7.4781598567 * t;
    fa[12] = 5.311886287 + 3.8133035638 * t;
    fa[13] = (0.024381750 + 0.00000538691 * t) * t;
    for slot in fa.iter_mut() {
        *slot = floored_mod(*slot, FRAME_TAU);
    }
    let ke1_f: Vec<f64> = IAU2000_KE1.iter().map(|&m| m as f64).collect();
    let a = series_dot(&ke1_f, &fa)?;
    let mut c_terms = IAU2000_SE1_0 * scalar_sin(a);
    c_terms += IAU2000_SE1_1 * scalar_cos(a);
    c_terms *= t;
    let ke0 = integer_matrix_vector(IAU2000_KE0_T.iter(), IAU2000_KE0_T.len(), 14, &fa)?;
    let mut sin_ke0 = Vec::with_capacity(ke0.len());
    let mut cos_ke0 = Vec::with_capacity(ke0.len());
    for &a in ke0.iter() {
        sin_ke0.push(scalar_sin(a));
        cos_ke0.push(scalar_cos(a));
    }
    c_terms += series_dot(&IAU2000_SE0_T0, &sin_ke0)?;
    c_terms += series_dot(&IAU2000_SE0_T1, &cos_ke0)?;
    Ok(c_terms * FRAME_ASEC2RAD)
}

/// Capitaine P03 dimensionless precession matrix at a finite TDB Julian date.
/// It rotates J2000 mean coordinates to the mean equator and equinox of date.
/// Non-finite dates return an error.
pub fn precession_matrix(jd_tdb: f64) -> Result<[[f64; 3]; 3], Error> {
    let jd_tdb = check_finite(jd_tdb)?;
    let delta = jd_tdb - FRAME_EPOCH_J2000;
    let t = delta / 36525.0;
    let value = -0.0000000951 * t;
    let value = value + 0.000132851;
    let value = value * t;
    let value = value - 0.00114045;
    let value = value * t;
    let value = value - 1.0790069;
    let value = value * t;
    let value = value + 5038.481507;
    let psia = value * t;
    let value = 0.0000003337 * t;
    let value = value - 0.000000467;
    let value = value * t;
    let value = value - 0.00772503;
    let value = value * t;
    let value = value + 0.0512623;
    let value = value * t;
    let value = value - 0.025754;
    let value = value * t;
    let omega = value + 84381.406;
    let value = -0.0000000560 * t;
    let value = value + 0.000170663;
    let value = value * t;
    let value = value - 0.00121197;
    let value = value * t;
    let value = value - 2.3814292;
    let value = value * t;
    let value = value + 10.556403;
    let chia = value * t;
    let eps = 84381.406 * FRAME_ASEC2RAD;
    let psia = psia * FRAME_ASEC2RAD;
    let omega = omega * FRAME_ASEC2RAD;
    let chia = chia * FRAME_ASEC2RAD;
    let (sa, ca) = eps.sin_cos();
    let (sb, cb) = (-psia).sin_cos();
    let (sc, cc) = (-omega).sin_cos();
    let (sd, cd) = chia.sin_cos();
    Ok([
        [
            cd * cb - sb * sd * cc,
            cd * sb * ca + sd * cc * cb * ca - sa * sd * sc,
            cd * sb * sa + sd * cc * cb * sa + ca * sd * sc,
        ],
        [
            -sd * cb - sb * cd * cc,
            -sd * sb * ca + cd * cc * cb * ca - sa * cd * sc,
            -sd * sb * sa + cd * cc * cb * sa + ca * cd * sc,
        ],
        [sb * sc, -sc * cb * ca - sa * cc, -sc * cb * sa + cc * ca],
    ])
}

/// Dimensionless nutation rotation from mean to true equator and equinox.
/// Mean obliquity, true obliquity and longitude nutation are radians.
/// Non-finite angles return an error.
pub fn nutation_matrix(
    mean_obliquity_rad: f64,
    true_obliquity_rad: f64,
    psi_rad: f64,
) -> Result<[[f64; 3]; 3], Error> {
    let mean_obliquity_rad = check_finite(mean_obliquity_rad)?;
    let true_obliquity_rad = check_finite(true_obliquity_rad)?;
    let psi_rad = check_finite(psi_rad)?;
    let (sobm, cobm) = mean_obliquity_rad.sin_cos();
    let (sobt, cobt) = true_obliquity_rad.sin_cos();
    let (spsi, cpsi) = psi_rad.sin_cos();
    Ok([
        [cpsi, -spsi * cobm, -spsi * sobm],
        [
            spsi * cobt,
            cpsi * cobm * cobt + sobm * sobt,
            cpsi * sobm * cobt - cobm * sobt,
        ],
        [
            spsi * sobt,
            cpsi * cobm * sobt - sobm * cobt,
            cpsi * sobm * sobt + cobm * cobt,
        ],
    ])
}

/// ICRS-to-J2000 dimensionless frame-bias matrix using IERS 2003 offsets.
/// Second-order diagonal terms preserve orthogonality at this angular scale.
pub fn frame_bias_matrix() -> [[f64; 3]; 3] {
    let xi0 = -0.0166170 * FRAME_ASEC2RAD;
    let eta0 = -0.0068192 * FRAME_ASEC2RAD;
    let da0 = -0.01460 * FRAME_ASEC2RAD;
    let yx = -da0;
    let zx = xi0;
    let xy = da0;
    let zy = eta0;
    let xz = -xi0;
    let yz = -eta0;
    let xx = 1.0 - 0.5 * (yx * yx + zx * zx);
    let yy = 1.0 - 0.5 * (yx * yx + zy * zy);
    let zz = 1.0 - 0.5 * (zy * zy + zx * zx);
    [[xx, xy, xz], [yx, yy, yz], [zx, zy, zz]]
}

/// Matrix-vector product in a dimensionless orthonormal frame.
/// Non-finite matrix entries or coordinates return an error.
pub fn matvec3(m: [[f64; 3]; 3], v: [f64; 3]) -> Result<[f64; 3], Error> {
    crate::positions::mat3_vec3(m, v)
}

/// Transpose of a 3x3 matrix (exact, no rounding).
pub fn transpose3(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [m[0][0], m[1][0], m[2][0]],
        [m[0][1], m[1][1], m[2][1]],
        [m[0][2], m[1][2], m[2][2]],
    ]
}

/// Passive x-axis rotation by `theta` radians. Non-finite inputs yield NaN.
pub fn rot_x_neg(theta: f64) -> [[f64; 3]; 3] {
    let c = scalar_cos(-theta);
    let s = scalar_sin(-theta);
    [[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]]
}

/// Row-major 3×3 matrix product. Non-finite entries return an error.
pub fn matmul3(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> Result<[[f64; 3]; 3], Error> {
    crate::positions::mat3_mul(a, b)
}

/// ICRS-to-equinox-of-date matrix, `N × P × B`.
/// Inputs are dimensionless nutation and precession matrices; non-finite
/// entries return an error. The bias uses IERS 2003 angular offsets.
pub fn equinox_matrix(
    nutation: [[f64; 3]; 3],
    precession: [[f64; 3]; 3],
) -> Result<[[f64; 3]; 3], Error> {
    matmul3(matmul3(nutation, precession)?, frame_bias_matrix())
}

/// Earth rotation angle in turns at a two-part UT1 Julian date, IAU 2000.
///
/// Both operands are days; they remain separate through the compensated ERA
/// law shared with the scalar-JD rotation entry. The second part need not be
/// confined to one day. The result is reduced to `[0,1)`; nonfinite operands
/// return an error. Finite overflow retains IEEE arithmetic behavior, and no
/// precession, nutation or Delta-T model is acquired here.
pub fn earth_rotation_angle_turns(jd_whole: f64, fraction_ut1: f64) -> Result<f64, Error> {
    let jd_whole = check_finite(jd_whole)?;
    let fraction_ut1 = check_finite(fraction_ut1)?;
    Ok(super::earth_rotation_turns_from_ut1_parts(
        jd_whole,
        fraction_ut1,
    ))
}

/// Greenwich mean sidereal time in hours, using the equinox method.
/// The TT Julian date and the UT1 fractional day must be finite.
/// The UT1 fraction shares the floored whole day of TT.
pub fn gmst_hours(jd_tt: f64, ut1_fraction: f64) -> Result<f64, Error> {
    let jd_tt = check_finite(jd_tt)?;
    let ut1_fraction = check_finite(ut1_fraction)?;
    let whole = jd_tt.floor();
    gmst_hours_two_part(whole, jd_tt - whole, ut1_fraction)
}

/// Mean sidereal time in hours with independent whole-TT and fraction inputs.
///
/// Both fractions are days; the UT1 fraction is supplied independently of
/// the TT fraction. Earth rotation and TDB precession-in-RA preserve this
/// representation instead of quantizing a short sample offset to a scalar JD.
/// Inputs must be finite. The returned angle is reduced to [0, 24) hours.
pub(crate) fn gmst_hours_two_part(
    whole: f64,
    tt_fraction: f64,
    ut1_fraction: f64,
) -> Result<f64, Error> {
    let whole = check_finite(whole)?;
    let tt_fraction = check_finite(tt_fraction)?;
    let ut1_fraction = check_finite(ut1_fraction)?;
    let theta = earth_rotation_angle_turns(whole, ut1_fraction)?;
    let tdb_fraction = tt_fraction + tdb_minus_tt_seconds(whole, tt_fraction)? / FRAME_DAY_S;

    let t = (whole - FRAME_EPOCH_J2000 + tdb_fraction) / 36525.0;
    let st = 0.014506
        + ((((-0.0000000368 * t - 0.000029956) * t - 0.00000044) * t + 1.3915817) * t
            + 4612.156534)
            * t;
    Ok(floored_mod(st / 54000.0 + theta * 24.0, 24.0))
}

/// Greenwich apparent sidereal time in hours, reduced to [0, 24).
/// Inputs are the TT Julian date, TT minus UT1 in seconds and the TDB Julian
/// date. The time-scale offset stays separate to preserve sub-day precision.
/// Non-finite inputs return an error.
pub fn gast_hours(jd_tt: f64, delta_t_seconds: f64, jd_tdb: f64) -> Result<f64, Error> {
    let jd_tt = check_finite(jd_tt)?;
    let whole = jd_tt.floor();
    gast_hours_two_part(whole, jd_tt - whole, delta_t_seconds, jd_tdb)
}

/// Apparent sidereal time in hours at a fraction-preserving TT observation.
///
/// `whole` and `tt_fraction` are days; TT-minus-UT1 is seconds and `jd_tdb`
/// is scalar TDB for mean obliquity. The rotation clock retains the TT fraction,
/// whereas nutation/complementary terms use the scalar TT observation, matching
/// their separate acquisition routes. Inputs and collapsed TT must be finite;
/// returns [0, 24) hours. No polar-motion or caller-model substitution occurs.
pub(crate) fn gast_hours_two_part(
    whole: f64,
    tt_fraction: f64,
    delta_t_seconds: f64,
    jd_tdb: f64,
) -> Result<f64, Error> {
    let whole = check_finite(whole)?;
    let tt_fraction = check_finite(tt_fraction)?;
    let jd_tt = check_finite(whole + tt_fraction)?;
    let delta_t_seconds = check_finite(delta_t_seconds)?;
    let jd_tdb = check_finite(jd_tdb)?;
    let ut1_fraction = tt_fraction - delta_t_seconds / FRAME_DAY_S;
    let gmst = gmst_hours_two_part(whole, tt_fraction, ut1_fraction)?;
    let (dpsi_u, _) = iau2000a_tenth_microarcsec(jd_tt)?;
    let dpsi = dpsi_u * FRAME_TENTH_USEC_2_RAD;
    let mean_obl = mean_obliquity_arcsec(jd_tdb)? * FRAME_ASEC2RAD;
    let c_terms = equinoxes_complimentary_terms_radians(jd_tt)?;
    let eq_eq = dpsi * mean_obl.cos() + c_terms;
    Ok(floored_mod(gmst + eq_eq / FRAME_TAU * 24.0, 24.0))
}

/// WGS84 geodetic position in ITRS, in astronomical units.
/// Inputs are east longitude and latitude in degrees and height in metres.
/// Non-finite inputs return an error.
pub fn itrs_position_au(lon_deg: f64, lat_deg: f64, height_m: f64) -> Result<[f64; 3], Error> {
    let lon_deg = check_finite(lon_deg)?;
    let lat_deg = check_finite(lat_deg)?;
    let height_m = check_finite(height_m)?;

    let lat = lat_deg / 360.0 * FRAME_TAU;
    let lon = lon_deg / 360.0 * FRAME_TAU;
    let radius_au = FRAME_WGS84_A_M / FRAME_AU_M;
    let elevation_au = height_m / FRAME_AU_M;
    let sinphi = scalar_sin(lat);
    let cosphi = scalar_cos(lat);
    let omf = (FRAME_WGS84_INV_F - 1.0) / FRAME_WGS84_INV_F;
    let omf2 = omf * omf;
    let c = 1.0 / (cosphi * cosphi + sinphi * sinphi * omf2).sqrt();
    let s = omf2 * c;
    let radius_xy = radius_au * c;
    let xy = (radius_xy + elevation_au) * cosphi;
    let x = xy * scalar_cos(lon);
    let y = xy * scalar_sin(lon);
    let radius_z = radius_au * s;
    let z = (radius_z + elevation_au) * sinphi;
    Ok([x, y, z])
}

/// ITRS velocity in AU/day from the nominal Earth spin, `ω × r`.
/// The position is in astronomical units. Non-finite inputs return an error.
pub fn itrs_velocity_au_per_day(itrs_au: [f64; 3]) -> Result<[f64; 3], Error> {
    for &v in itrs_au.iter() {
        check_finite(v)?;
    }
    let k = FRAME_ANGVEL_RAD_S * FRAME_DAY_S;
    Ok([k * (-itrs_au[1]), k * itrs_au[0], k * (0.0 * itrs_au[2])])
}

/// Celestial-to-terrestrial rotation without polar motion.
/// Inputs are the dimensionless equinox matrix and GAST in hours.
/// Non-finite entries or angles return an error.
pub fn itrs_rotation_matrix(equinox: [[f64; 3]; 3], gast_h: f64) -> Result<[[f64; 3]; 3], Error> {
    let gast_h = check_finite(gast_h)?;
    for row in equinox.iter() {
        for &v in row.iter() {
            check_finite(v)?;
        }
    }
    let spin = rot_z(-gast_h * FRAME_TAU / 24.0);
    matmul3(spin, equinox)
}

/// Geocentric GCRS site position and velocity in AU and AU/day.
/// Inputs are the dimensionless celestial-to-terrestrial matrix and an ITRS
/// position in AU. The transpose supplies the inverse rotation.
/// Non-finite coordinates or matrix entries return an error.
pub fn itrs_to_gcrs(
    rotation: [[f64; 3]; 3],
    itrs_au: [f64; 3],
) -> Result<([f64; 3], [f64; 3]), Error> {
    let rt = transpose3(rotation);
    let v = itrs_velocity_au_per_day(itrs_au)?;
    Ok((matvec3(rt, itrs_au)?, matvec3(rt, v)?))
}

/// True-ecliptic-of-date dimensionless rotation.
/// Inputs are the equinox matrix, mean obliquity and obliquity nutation in
/// radians. Non-finite entries or angles return an error.
pub fn ecliptic_date_matrix(
    equinox: [[f64; 3]; 3],
    mean_obliquity_rad: f64,
    deps_rad: f64,
) -> Result<[[f64; 3]; 3], Error> {
    let mean_obliquity_rad = check_finite(mean_obliquity_rad)?;
    let deps_rad = check_finite(deps_rad)?;
    let rx = rot_x_neg(mean_obliquity_rad + deps_rad);
    matmul3(rx, equinox)
}

/// Spherical coordinates in a dimensionless reference frame.
/// Returns longitude in [0, 360) degrees, latitude in degrees and distance in
/// AU. Coordinates must be finite; the zero vector has zero distance.
pub fn frame_latlon_deg(frame: [[f64; 3]; 3], xyz_au: [f64; 3]) -> Result<[f64; 3], Error> {
    for row in frame.iter() {
        for &v in row.iter() {
            check_finite(v)?;
        }
    }
    for &v in xyz_au.iter() {
        check_finite(v)?;
    }
    let v = matvec3(frame, xyz_au)?;
    let r = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    let theta = v[2].atan2((v[0] * v[0] + v[1] * v[1]).sqrt());
    let phi = floored_mod(v[1].atan2(v[0]), FRAME_TAU);
    let lon = floored_mod(phi * 360.0 / FRAME_TAU, 360.0);
    let lat = theta * 360.0 / FRAME_TAU;
    Ok([lon, lat, r])
}

/// Barycentric observer position, velocity and geocentric site offset.
/// Units are AU, AU/day and AU, respectively.
pub type BarycentricSiteState = ([f64; 3], [f64; 3], [f64; 3]);

/// Barycentric observer position and velocity with a WGS84 topocentric offset.
/// Earth position and velocity are AU and AU/day. Longitude/latitude are
/// degrees, height is metres, epochs are Julian days and Delta-T is seconds.
/// Nutation, precession and Earth spin supply the rotation without polar
/// motion. Non-finite inputs return an error.
#[allow(clippy::too_many_arguments)]
pub fn barycentric_site_state(
    earth_pos: [f64; 3],
    earth_vel: [f64; 3],
    lon_deg: f64,
    lat_deg: f64,
    height_m: f64,
    jd_tt: f64,
    delta_t_seconds: f64,
    jd_tdb: f64,
) -> Result<BarycentricSiteState, Error> {
    for &v in earth_pos.iter().chain(earth_vel.iter()) {
        check_finite(v)?;
    }
    let (dpsi_u, deps_u) = iau2000a_tenth_microarcsec(jd_tt)?;
    let dpsi = dpsi_u * FRAME_TENTH_USEC_2_RAD;
    let deps = deps_u * FRAME_TENTH_USEC_2_RAD;
    let mean_obl = mean_obliquity_arcsec(jd_tdb)? * FRAME_ASEC2RAD;
    let nut = nutation_matrix(mean_obl, mean_obl + deps, dpsi)?;
    let prec = precession_matrix(jd_tdb)?;
    let equinox = equinox_matrix(nut, prec)?;
    let gast = gast_hours(jd_tt, delta_t_seconds, jd_tdb)?;
    let rotation = itrs_rotation_matrix(equinox, gast)?;
    let itrs = itrs_position_au(lon_deg, lat_deg, height_m)?;
    let (off_pos, off_vel) = itrs_to_gcrs(rotation, itrs)?;
    let pos = [
        earth_pos[0] + off_pos[0],
        earth_pos[1] + off_pos[1],
        earth_pos[2] + off_pos[2],
    ];
    let vel = [
        earth_vel[0] + off_vel[0],
        earth_vel[1] + off_vel[1],
        earth_vel[2] + off_vel[2],
    ];
    Ok((pos, vel, off_pos))
}

/// Celestial-to-terrestrial rotation at a TT Julian date, without polar motion.
/// The bundled TT-minus-UT1 table supplies the independent rotation clock;
/// `jd_ut1` identifies the corresponding epoch and must be finite. For an
/// explicit independent UT1 offset, use [`itrs_rotation_at_parts`].
/// All non-finite inputs return an error.
pub fn itrs_rotation_at_tt_ut1(jd_tt: f64, jd_ut1: f64) -> Result<[[f64; 3]; 3], Error> {
    let jd_tt = check_finite(jd_tt)?;
    let _jd_ut1 = check_finite(jd_ut1)?;
    let whole = jd_tt.floor();
    let tt_fraction = jd_tt - whole;
    let delta_t_seconds = crate::smh2016::bundled_deltat_seconds(jd_tt)?;
    let ut1_fraction = tt_fraction - delta_t_seconds / FRAME_DAY_S;
    rotation_from_whole_fractions(whole, tt_fraction, ut1_fraction, jd_tt)
}

/// Celestial-to-terrestrial rotation from independent whole and fractional days.
/// `tt_fraction` and `ut1_fraction` share `whole`. Scalar `jd_tt` drives
/// nutation; the split TT epoch drives TDB precession and sidereal time.
/// Keeping UT1 separate avoids cancellation at Julian-date magnitude.
/// Polar motion is excluded. Non-finite inputs return an error.
pub fn itrs_rotation_at_parts(
    whole: f64,
    tt_fraction: f64,
    ut1_fraction: f64,
    jd_tt: f64,
) -> Result<[[f64; 3]; 3], Error> {
    let whole = check_finite(whole)?;
    let tt_fraction = check_finite(tt_fraction)?;
    let ut1_fraction = check_finite(ut1_fraction)?;
    let jd_tt = check_finite(jd_tt)?;
    rotation_from_whole_fractions(whole, tt_fraction, ut1_fraction, jd_tt)
}

/// Shared equinox-chain tail over an explicit day split: TDB, IAU 2000A
/// nutation at scalar TT, P03 precession/mean-obliquity at TDB, GAST from
/// the given TT and UT1 fractions, and `Rz(−GAST) @ M`. The GAST body
/// mirrors [`gast_hours_two_part`] except the UT1 fraction arrives
/// explicitly instead of being recovered from the TT fraction minus the
/// table offset. No polar motion. All inputs must be finite.
fn rotation_from_whole_fractions(
    whole: f64,
    tt_fraction: f64,
    ut1_fraction: f64,
    jd_tt: f64,
) -> Result<[[f64; 3]; 3], Error> {
    let (equinox, gast) = equinox_and_gast_at_parts(whole, tt_fraction, ut1_fraction, jd_tt)?;
    itrs_rotation_matrix(equinox, gast)
}

/// True equator/equinox-of-date matrix and GAST in hours.
/// TT and UT1 fractions share one whole Julian day. Nutation uses scalar TT;
/// precession and mean obliquity use the corresponding TDB epoch.
/// No polar motion is applied. Non-finite inputs return an error.
pub(crate) fn equinox_and_gast_at_parts(
    whole: f64,
    tt_fraction: f64,
    ut1_fraction: f64,
    jd_tt: f64,
) -> Result<([[f64; 3]; 3], f64), Error> {
    let whole = check_finite(whole)?;
    let tt_fraction = check_finite(tt_fraction)?;
    let ut1_fraction = check_finite(ut1_fraction)?;
    let jd_tt = check_finite(jd_tt)?;
    // Add short TDB and TT intervals before collapsing the full epoch.
    let jd_tdb = check_finite(
        whole + (tt_fraction + tdb_minus_tt_seconds(whole, tt_fraction)? / FRAME_DAY_S),
    )?;
    let (dpsi_u, deps_u) = iau2000a_tenth_microarcsec(jd_tt)?;
    let dpsi = dpsi_u * FRAME_TENTH_USEC_2_RAD;
    let deps = deps_u * FRAME_TENTH_USEC_2_RAD;
    let mean_obl = mean_obliquity_arcsec(jd_tdb)? * FRAME_ASEC2RAD;
    let nut = nutation_matrix(mean_obl, mean_obl + deps, dpsi)?;
    let prec = precession_matrix(jd_tdb)?;
    let equinox = equinox_matrix(nut, prec)?;
    let gmst = gmst_hours_two_part(whole, tt_fraction, ut1_fraction)?;
    let c_terms = equinoxes_complimentary_terms_radians(jd_tt)?;
    let eq_eq = dpsi * mean_obl.cos() + c_terms;
    let gast = floored_mod(gmst + eq_eq / FRAME_TAU * 24.0, 24.0);
    Ok((equinox, gast))
}

/// Independent TT and UT1 fractional days for one observation epoch.
/// The whole day is shared so sub-day offsets never require subtraction
/// between two full-magnitude Julian dates.
#[derive(Debug, Clone, Copy)]
pub struct VectorRotationEpochs {
    /// Shared whole Julian day.
    pub whole: f64,
    /// TT fraction in days relative to the shared whole day.
    pub tt_fraction: f64,
    /// UT1 fraction in days relative to the shared whole day.
    pub ut1_fraction: f64,
}

impl VectorRotationEpochs {
    /// TT/UT1 fractions from a finite UT1 Julian date using bundled Delta-T.
    /// Returns both the split epoch and the scalar TT Julian date.
    pub fn from_ut1(jd_ut1: f64) -> Result<(Self, f64), Error> {
        let (whole, tt_fraction) = crate::smh2016::timescale_parts_from_ut1(jd_ut1)?;
        Ok((
            Self {
                whole,
                tt_fraction,
                ut1_fraction: 0.0,
            },
            whole + tt_fraction,
        ))
    }

    /// TT/UT1 fractions from a finite TT Julian date using bundled Delta-T.
    /// Non-finite input dates return an error.
    pub fn from_tt(jd_tt: f64) -> Result<Self, Error> {
        let jd_tt = check_finite(jd_tt)?;
        let whole = jd_tt.floor();
        let tt_fraction = jd_tt - whole;
        let delta_t_seconds = crate::smh2016::bundled_deltat_seconds(jd_tt)?;
        Ok(Self {
            whole,
            tt_fraction,
            ut1_fraction: tt_fraction - delta_t_seconds / FRAME_DAY_S,
        })
    }
}

/// Light-time-corrected TT Julian date from a two-part TT observation epoch.
/// Converts the fraction to TDB, subtracts the light time in days, then
/// converts back to TT before collapsing the result. The split arithmetic
/// preserves precision for light times much smaller than one Julian day.
/// Non-finite epochs or intervals return an error.
pub fn split_retarded_epoch_tt(
    whole: f64,
    tt_fraction: f64,
    light_time_days: f64,
) -> Result<f64, Error> {
    let whole = check_finite(whole)?;
    let tt_fraction = check_finite(tt_fraction)?;
    let light_time_days = check_finite(light_time_days)?;
    let tdb_fraction = crate::time::tt_to_tdb_fraction(whole, tt_fraction)?;
    let retarded = check_finite(tdb_fraction - light_time_days)?;
    let back = crate::time::tdb_minus_tt_seconds(whole, retarded)? / FRAME_DAY_S;
    check_finite(whole + (retarded - back))
}

/// Active z-axis rotation by `theta` radians. Non-finite inputs yield NaN.
pub fn rot_z(theta: f64) -> [[f64; 3]; 3] {
    let c = scalar_cos(theta);
    let s = scalar_sin(theta);
    [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_rotation(matrix: [[f64; 3]; 3]) {
        let identity = matmul3(matrix, transpose3(matrix)).unwrap();
        for (r, row) in identity.iter().enumerate() {
            for (c, &value) in row.iter().enumerate() {
                let expected = if r == c { 1.0 } else { 0.0 };
                assert!((value - expected).abs() < 2e-13, "[{r}][{c}]={value}");
            }
        }
        let vector = [0.3, -0.4, 0.5];
        let rotated = matvec3(matrix, vector).unwrap();
        let inverse = matvec3(transpose3(matrix), rotated).unwrap();
        for (actual, expected) in inverse.into_iter().zip(vector) {
            assert!((actual - expected).abs() < 2e-13);
        }
    }

    #[test]
    fn precession_and_nutation_are_rotation_matrices() {
        for jd_tt in [
            1_500_000.5,
            2_378_496.5,
            FRAME_EPOCH_J2000,
            2_460_310.5,
            3_000_000.5,
        ] {
            let tdb = crate::time::tt_to_tdb_jd(jd_tt).unwrap();
            let epsilon = mean_obliquity_arcsec(tdb).unwrap() * FRAME_ASEC2RAD;
            let (dpsi, deps) = iau2000a_tenth_microarcsec(jd_tt).unwrap();
            let nutation = nutation_matrix(
                epsilon,
                epsilon + deps * FRAME_TENTH_USEC_2_RAD,
                dpsi * FRAME_TENTH_USEC_2_RAD,
            )
            .unwrap();
            let precession = precession_matrix(tdb).unwrap();
            assert_rotation(precession);
            assert_rotation(nutation);
            let equinox = equinox_matrix(nutation, precession).unwrap();
            assert_rotation(equinox);
            assert_rotation(
                ecliptic_date_matrix(equinox, epsilon, deps * FRAME_TENTH_USEC_2_RAD).unwrap(),
            );
        }
    }

    #[test]
    fn frame_origins_and_zero_nutation_have_their_defining_values() {
        let epsilon = 84_381.406 * FRAME_ASEC2RAD;
        assert!((mean_obliquity_arcsec(FRAME_EPOCH_J2000).unwrap() - 84_381.406).abs() < 1e-10);
        for matrix in [
            precession_matrix(FRAME_EPOCH_J2000).unwrap(),
            nutation_matrix(epsilon, epsilon, 0.0).unwrap(),
        ] {
            for (r, row) in matrix.iter().enumerate() {
                for (c, &value) in row.iter().enumerate() {
                    let expected = if r == c { 1.0 } else { 0.0 };
                    assert!((value - expected).abs() < 4.0 * f64::EPSILON);
                }
            }
        }
    }

    #[test]
    fn wgs84_sites_and_earth_spin_follow_geometric_definitions() {
        let a = FRAME_WGS84_A_M / FRAME_AU_M;
        let b = a * (1.0 - 1.0 / FRAME_WGS84_INV_F);
        let equator = itrs_position_au(0.0, 0.0, 0.0).unwrap();
        let pole = itrs_position_au(0.0, 90.0, 0.0).unwrap();
        assert!((equator[0] - a).abs() < 1e-19);
        assert!(equator[1].abs() < 1e-19 && equator[2].abs() < 1e-19);
        assert!(pole[0].abs() < 1e-19 && pole[1].abs() < 1e-19);
        assert!((pole[2] - b).abs() < 1e-19);
        let speed = itrs_velocity_au_per_day(equator).unwrap();
        let expected = a * FRAME_ANGVEL_RAD_S * FRAME_DAY_S;
        assert!(speed[0].abs() < 1e-19 && speed[2].abs() < 1e-19);
        assert!((speed[1] - expected).abs() < 1e-19);
        let rotation = itrs_rotation_matrix(frame_bias_matrix(), 7.0).unwrap();
        let (gcrs, velocity) = itrs_to_gcrs(rotation, equator).unwrap();
        assert!((series_dot(&gcrs, &gcrs).unwrap().sqrt() - a).abs() < 1e-18);
        assert!(series_dot(&gcrs, &velocity).unwrap().abs() < 1e-20);
    }

    #[test]
    fn split_clocks_keep_sub_day_offsets_and_rotation_periodicity() {
        for jd_ut1 in [2_378_496.5, FRAME_EPOCH_J2000, 2_460_310.5] {
            let (epochs, tt) = VectorRotationEpochs::from_ut1(jd_ut1).unwrap();
            assert_eq!(epochs.whole, jd_ut1);
            assert_eq!(epochs.ut1_fraction, 0.0);
            assert!(
                (epochs.tt_fraction * FRAME_DAY_S
                    - crate::smh2016::bundled_deltat_seconds(tt).unwrap())
                .abs()
                    < 1e-5
            );
            let frame =
                itrs_rotation_at_parts(epochs.whole, epochs.tt_fraction, epochs.ut1_fraction, tt)
                    .unwrap();
            assert_rotation(frame);
            let plus_turn = itrs_rotation_matrix(frame, 24.0).unwrap();
            for (actual, expected) in plus_turn.iter().flatten().zip(frame.iter().flatten()) {
                assert!((actual - expected).abs() < 1e-15);
            }
        }
    }

    #[test]
    fn zero_light_time_and_spherical_projection_preserve_distance() {
        let whole = FRAME_EPOCH_J2000;
        let fraction = 0.25;
        let unretarded = split_retarded_epoch_tt(whole, fraction, 0.0).unwrap();
        assert!((unretarded - (whole + fraction)).abs() <= f64::EPSILON * whole);
        let light_time = 0.001;
        let retarded = split_retarded_epoch_tt(whole, fraction, light_time).unwrap();
        assert!(((unretarded - retarded) - light_time).abs() < 1e-9);
        let xyz = [0.3, -0.4, 0.5];
        let place = frame_latlon_deg(rot_z(0.8), xyz).unwrap();
        assert!((place[2] - 0.5f64.sqrt()).abs() < 1e-15);
        assert!(place[0] >= 0.0 && place[0] < 360.0);
        assert!(place[1] >= -90.0 && place[1] <= 90.0);
    }

    #[test]
    fn scalar_series_products_and_invalid_input_domains() {
        let weights = [0.2, -0.1, 0.4];
        let matrix = [[1, 2, 3], [4, 5, 6]];
        let product = integer_matrix_vector(matrix, 2, 3, &weights).unwrap();
        let flat = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        assert_eq!(product, matrix_vector(&flat, 2, 3, &weights).unwrap());
        assert!((product[0] - 1.2).abs() < 1e-15);
        assert!((product[1] - 2.7).abs() < 1e-15);
        assert!(series_dot(&[f64::NAN], &[1.0]).is_err());
        assert!(precession_matrix(f64::INFINITY).is_err());
        assert!(iau2000a_tenth_microarcsec(f64::NAN).is_err());
        assert!(equinoxes_complimentary_terms_radians(f64::INFINITY).is_err());
        assert!(itrs_rotation_at_parts(FRAME_EPOCH_J2000, 0.0, 0.0, f64::NAN).is_err());
        assert!(split_retarded_epoch_tt(FRAME_EPOCH_J2000, 0.0, f64::NAN).is_err());
    }
}

#[cfg(test)]
mod exact_nutation_cache_tests {
    use super::*;

    #[test]
    fn cached_full_series_preserves_both_angle_words() {
        for epoch in [1000000.25, 2415020.5, 2451545.0, 2488069.5, 5000000.75] {
            let expected = iau2000a_uncached(epoch).unwrap();
            for _ in 0..3 {
                let actual = iau2000a_tenth_microarcsec(epoch).unwrap();
                assert_eq!(
                    [actual.0.to_bits(), actual.1.to_bits()],
                    [expected.0.to_bits(), expected.1.to_bits()]
                );
            }
        }
        assert!(iau2000a_tenth_microarcsec(f64::INFINITY).is_err());
    }
}
