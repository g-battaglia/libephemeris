//! Constant sidereal frame construction retaining physical arithmetic residuals.
//!
//! Scalar positions keep their existing public transformation. Daily rates
//! carry every original model term through the same fixed projection before
//! the final Jacobian; no matrix entry is promoted to an exact physical datum.
use super::{derivative::AngularWord as Word, *};
use crate::constants::{SIDBIT_ECL_T0, SIDM_GALALIGN_MARDYKS, SIDM_J2000};
use crate::sidereal;

pub(crate) type WordMatrix = [[Word; 3]; 3];

/// Already acquired operands of the mean epoch-pair zero direction.
#[derive(Clone, Copy)]
pub(crate) struct EpochPairZero {
    pub(crate) tt: f64,
    pub(crate) defining_tt: f64,
    pub(crate) defining_degrees: f64,
    pub(crate) ecliptic_date: bool,
}

fn identity() -> WordMatrix {
    std::array::from_fn(|i| std::array::from_fn(|j| Word::from(if i == j { 1.0 } else { 0.0 })))
}
fn dot(a: [Word; 3], b: [Word; 3]) -> Word {
    a.into_iter()
        .zip(b)
        .fold(Word::from(0.0), |sum, (a, b)| sum + a * b)
}
fn cross(a: [Word; 3], b: [Word; 3]) -> [Word; 3] {
    std::array::from_fn(|i| a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3])
}
fn unit(v: [Word; 3]) -> Result<[Word; 3]> {
    let length = dot(v, v).sqrt();
    if length.rounded() == 0.0 || !length.rounded().is_finite() {
        return Err(Error::MathDomain {
            value: length.rounded(),
        });
    }
    Ok(v.map(|v| v / length))
}
fn transpose(m: WordMatrix) -> WordMatrix {
    std::array::from_fn(|i| std::array::from_fn(|j| m[j][i]))
}
fn multiply(a: WordMatrix, b: WordMatrix) -> WordMatrix {
    std::array::from_fn(|i| std::array::from_fn(|j| dot(a[i], [b[0][j], b[1][j], b[2][j]])))
}
fn radians(degrees: f64) -> Word {
    // The existing verified two-word mathematical degree conversion.
    Word::from(degrees)
        / Word {
            hi: 57.29577951308232,
            lo: -1.9878495670576283e-15,
        }
}
fn rx(angle: Word) -> Result<WordMatrix> {
    let (s, c) = super::frame_velocity::frame_sin_cos_words(angle)?;
    let z = Word::from(0.0);
    let o = Word::from(1.0);
    Ok([[o, z, z], [z, c, s], [z, -s, c]])
}
fn rz(angle: Word) -> Result<WordMatrix> {
    let (s, c) = super::frame_velocity::frame_sin_cos_words(angle)?;
    Ok(rz_components(s, c))
}
fn rz_components(s: Word, c: Word) -> WordMatrix {
    let z = Word::from(0.0);
    let o = Word::from(1.0);
    [[c, s, z], [-s, c, z], [z, z, o]]
}

/// Complete published pole pairs, sharing the compensated series evaluator.
fn poles(t: Word, eps: Word) -> Result<([Word; 3], [Word; 3])> {
    let (p, q) = super::frame_velocity::vondrak_pair_words(
        t,
        &VONDRAK_PQPER_PERIODS,
        &VONDRAK_PQPER_P_COS,
        &VONDRAK_PQPER_P_SIN,
        &VONDRAK_PQPER_Q_COS,
        &VONDRAK_PQPER_Q_SIN,
        &VONDRAK_PQPOL,
    )?;
    let (x, y) = super::frame_velocity::vondrak_pair_words(
        t,
        &VONDRAK_XYPER_PERIODS,
        &VONDRAK_XYPER_X_COS,
        &VONDRAK_XYPER_X_SIN,
        &VONDRAK_XYPER_Y_COS,
        &VONDRAK_XYPER_Y_SIN,
        &VONDRAK_XYPOL,
    )?;
    let z = Word::from(0.0);
    let o = Word::from(1.0);
    let ordinate = |a: Word, b: Word| {
        let square = o - a * a - b * b;
        if square.rounded() > 0.0 {
            square.sqrt()
        } else {
            z
        }
    };
    let (s, c) = super::frame_velocity::frame_sin_cos_words(eps)?;
    let w = ordinate(p, q);
    Ok(([x, y, ordinate(x, y)], [p, -q * c - w * s, -q * s + w * c]))
}
fn epoch_poles(tt: f64) -> Result<([Word; 3], [Word; 3])> {
    // The Julian-epoch model map is evaluated before intermediate rounding.
    let epoch = Word::from(2000.0) + (Word::from(tt) - Word::from(2451545.0)) / Word::from(365.25);
    let t = (epoch - Word::from(2000.0)) / Word::from(100.0);
    poles(t, Word::from(84381.406) * Word::from(crate::erfa::DAS2R))
}
fn ecliptic_basis(tt: f64) -> Result<WordMatrix> {
    let (equator, ecliptic) = epoch_poles(tt)?;
    let x = unit(cross(ecliptic, equator))?;
    Ok([x, cross(ecliptic, x), ecliptic])
}
fn biased_equator(tt: f64) -> Result<WordMatrix> {
    let (equator, ecliptic) = epoch_poles(tt)?;
    let x = unit(cross(equator, ecliptic))?;
    let basis = [x, cross(equator, x), equator];
    let scale = Word::from(crate::erfa::DAS2R);
    let dx = Word::from(-0.016617) * scale;
    let de = Word::from(-0.0068192) * scale;
    let dr = Word::from(-0.0146) * scale;
    Ok(basis.map(|[x, y, z]| {
        [
            x - y * dr + z * dx,
            x * dr + y + z * de,
            -x * dx - y * de + z,
        ]
    }))
}
/// Fixed mean-equator transform; the defining J2000 identity is explicit.
pub(crate) fn equatorial_epoch(tt: f64) -> Result<WordMatrix> {
    if tt == sidereal::J2000_TT {
        return Ok(identity());
    }
    Ok(multiply(
        biased_equator(tt)?,
        transpose(biased_equator(sidereal::J2000_TT)?),
    ))
}
fn iau_obliquity(tt: f64) -> Word {
    let t = (Word::from(tt) - Word::from(2451545.0)) / Word::from(36525.0);
    [
        84381.406,
        -46.836769,
        -0.0001831,
        0.00200340,
        -0.000000576,
        -0.0000000434,
    ]
    .into_iter()
    .rev()
    .fold(Word::from(0.0), |sum, c| sum * t + Word::from(c))
        * Word::from(crate::erfa::DAS2R)
}
fn pole_obliquity_rotation(tt: f64) -> Result<WordMatrix> {
    let t = (Word::from(tt) - Word::from(2451545.0)) / Word::from(36525.0);
    let (equator, ecliptic) = poles(t, Word::from(VONDRAK_EPS0_RAD))?;
    // Preserve the scalar pole-angle acos clamp, including extrapolation.
    let one = Word::from(1.0);
    let raw = dot(equator, ecliptic);
    let c = if (raw - one).rounded() > 0.0 {
        one
    } else if (raw + one).rounded() < 0.0 {
        -one
    } else {
        raw
    };
    let squared = one - c * c;
    let s = if squared.rounded() > 0.0 {
        squared.sqrt()
    } else {
        Word::from(0.0)
    };
    let z = Word::from(0.0);
    let o = Word::from(1.0);
    Ok([[o, z, z], [z, c, s], [z, -s, c]])
}
/// Fixed-mode ecliptic/equatorial transform with its original obliquity model.
pub(crate) fn fixed_epoch(mode: i32, equatorial: bool) -> Result<WordMatrix> {
    let tt = sidereal::fixed_epoch_t0_jd(mode).ok_or(Error::AyanamshaAnchorUnserved { mode })?;
    let eq = equatorial_epoch(tt)?;
    if equatorial || mode == SIDM_J2000 {
        return Ok(eq);
    }
    let mut m = multiply(
        multiply(rx(iau_obliquity(tt))?, eq),
        rx(-iau_obliquity(sidereal::J2000_TT))?,
    );
    if mode == SIDM_GALALIGN_MARDYKS {
        m = multiply(rz(radians(sidereal::MARDYKS_DEFINING.0))?, m);
    }
    Ok(m)
}
fn zero_rotation(zero: EpochPairZero) -> Result<WordMatrix> {
    let defining = ecliptic_basis(zero.defining_tt)?;
    let current = ecliptic_basis(zero.tt)?;
    let (s0, c0) = super::frame_velocity::frame_sin_cos_words(radians(zero.defining_degrees))?;
    let (s, c) = if zero.ecliptic_date {
        let direction = std::array::from_fn(|i| defining[0][i] * c0 + defining[1][i] * s0);
        let x = dot(current[0], direction);
        let y = dot(current[1], direction);
        let length = (x * x + y * y).sqrt();
        (y / length, x / length)
    } else {
        let x = dot(defining[0], current[0]);
        let y = dot(defining[1], current[0]);
        let length = (x * x + y * y).sqrt();
        ((s0 * x - c0 * y) / length, (c0 * x + s0 * y) / length)
    };
    Ok(rz_components(s, c))
}
/// Ecliptic-of-T0 or invariable-plane projection, with the actual epoch-pair zero.
pub(crate) fn ecliptic_projection(bits: i32, tt: f64, zero: EpochPairZero) -> Result<WordMatrix> {
    let rotation = zero_rotation(zero)?;
    if bits & SIDBIT_ECL_T0 != 0 {
        let target = pole_obliquity_rotation(tt)?;
        let origin = transpose(pole_obliquity_rotation(sidereal::J2000_TT)?);
        return Ok(multiply(
            rotation,
            multiply(multiply(target, equatorial_epoch(tt)?), origin),
        ));
    }
    let plane = multiply(rx(radians(1.578690))?, rz(radians(107.582222))?);
    let direction = [rotation[0][0], rotation[0][1], Word::from(0.0)];
    let x = dot(plane[0], direction);
    let y = dot(plane[1], direction);
    let length = (x * x + y * y).sqrt();
    Ok(multiply(rz_components(y / length, x / length), plane))
}
/// Reconstruct a supplied spherical tuple without another binary64 reduction.
/// Undefined/zero-distance sentinels retain their existing scalar policy.
pub(crate) fn spherical_rates(xx: [f64; 6], matrix: WordMatrix) -> Result<Option<[f64; 3]>> {
    if xx[2] == 0.0 || xx.iter().any(|v| !v.is_finite()) {
        return Ok(None);
    }
    let (sl, cl) = super::frame_velocity::frame_sin_cos_words(radians(xx[0]))?;
    let (sb, cb) = super::frame_velocity::frame_sin_cos_words(radians(xx[1]))?;
    let r = Word::from(xx[2]);
    let dr = Word::from(xx[5]);
    let dl = radians(xx[3]);
    let db = radians(xx[4]);
    let p = [r * cb * cl, r * cb * sl, r * sb];
    let v = [
        dr * cb * cl - r * sb * cl * db - r * cb * sl * dl,
        dr * cb * sl - r * sb * sl * db + r * cb * cl * dl,
        dr * sb + r * cb * db,
    ];
    let z = Word::from(0.0);
    let rates = super::derivative::expanded_output_rates(matrix, [[z; 3]; 3], p, v, false)?;
    Ok(Some([rates.0, rates.1, rates.2]))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn complete_fixed_models_match_independent_mp120_matrices() {
        {
            let actual = fixed_epoch(crate::constants::SIDM_B1950, false).unwrap();
            let expected = [
                [
                    [0.999925716738747, 4.990698843365445e-17],
                    [0.012188555133020204, 8.36167246317226e-19],
                    [1.1326762210743623e-05, -5.839859791431151e-22],
                ],
                [
                    [-0.012188556339829055, 1.1004619768742843e-19],
                    [0.9999257103574211, 2.603002555768885e-17],
                    [0.00011340441716608828, 6.753076484468633e-22],
                ],
                [
                    [-9.943684756496124e-06, 7.257685466866669e-23],
                    [-0.00011353405000220733, -4.649659852884616e-21],
                    [0.9999999935055828, 2.8826120618473545e-17],
                ],
            ];
            for (a, e) in actual
                .into_iter()
                .flatten()
                .zip(expected.into_iter().flatten())
            {
                let expected = Word { hi: e[0], lo: e[1] };
                assert!(
                    (a - expected).rounded().abs() < 1e-27,
                    "{}",
                    (a - expected).rounded()
                );
            }
        }
        {
            let actual = equatorial_epoch(2435553.5).unwrap();
            let expected = [
                [
                    [0.9999430415203057, -4.957104101728982e-17],
                    [0.009788421053728171, 6.82667040213578e-19],
                    [0.004254471579225215, 2.110289224727921e-19],
                ],
                [
                    [-0.009788420824139785, 1.0333562714084351e-20],
                    [0.9999520920433835, 2.0047229017933325e-17],
                    [-2.087683490847244e-05, -8.023337214181738e-22],
                ],
                [
                    [-0.004254472107431036, 3.760765170362198e-19],
                    [-2.0768912417530773e-05, 6.135357160773607e-22],
                    [0.9999909494769212, 2.1121987997833654e-17],
                ],
            ];
            for (a, e) in actual
                .into_iter()
                .flatten()
                .zip(expected.into_iter().flatten())
            {
                let expected = Word { hi: e[0], lo: e[1] };
                assert!(
                    (a - expected).rounded().abs() < 1e-27,
                    "{}",
                    (a - expected).rounded()
                );
            }
        }
        {
            let actual = ecliptic_projection(
                256,
                2435553.5,
                EpochPairZero {
                    tt: 2435553.5003636805,
                    defining_tt: 2451545.0,
                    defining_degrees: 23.85709166666667,
                    ecliptic_date: false,
                },
            )
            .unwrap();
            let expected = [
                [
                    [0.9145571054950318, 1.1914643189182062e-17],
                    [0.4044567943130193, -1.2306044430346497e-17],
                    [4.819346701161677e-05, -2.6015598678233604e-21],
                ],
                [
                    [-0.40445679709001253, -7.739512365709836e-18],
                    [0.9145571013620617, 5.326774795525417e-17],
                    [8.738371896395308e-05, -5.801818961286622e-21],
                ],
                [
                    [-8.732738645247633e-06, 6.858288185632572e-22],
                    [-9.940957639795937e-05, -2.230932454819329e-22],
                    [0.9999999950207492, -1.3048381000829044e-17],
                ],
            ];
            for (a, e) in actual
                .into_iter()
                .flatten()
                .zip(expected.into_iter().flatten())
            {
                let expected = Word { hi: e[0], lo: e[1] };
                assert!(
                    (a - expected).rounded().abs() < 1e-27,
                    "{}",
                    (a - expected).rounded()
                );
            }
        }
        {
            let actual = ecliptic_projection(
                512,
                2435553.5,
                EpochPairZero {
                    tt: 2451545.000738761,
                    defining_tt: 2451545.0,
                    defining_degrees: 23.85709166666667,
                    ecliptic_date: false,
                },
            )
            .unwrap();
            let expected = [
                [
                    [0.9141807518581585, 2.0043408166611745e-17],
                    [0.4043805539673501, 8.95086667924722e-18],
                    [-0.02738467646533409, -4.9251322602861285e-19],
                ],
                [
                    [-0.40445496149878646, 9.97695878062951e-18],
                    [0.9145529565152579, -3.1366035451398657e-17],
                    [0.003012282891095285, 1.046090774601915e-19],
                ],
                [
                    [0.02626284544879258, -1.6523511325078942e-18],
                    [0.008322097227252467, 6.951940796448845e-19],
                    [0.9996204307869427, 7.976196572545898e-18],
                ],
            ];
            for (a, e) in actual
                .into_iter()
                .flatten()
                .zip(expected.into_iter().flatten())
            {
                let expected = Word { hi: e[0], lo: e[1] };
                assert!(
                    (a - expected).rounded().abs() < 1e-27,
                    "{}",
                    (a - expected).rounded()
                );
            }
        }
    }
}
