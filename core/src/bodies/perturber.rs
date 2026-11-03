//! Secular elements of Jupiter, Saturn, Uranus and Neptune.
//!
//! The coefficient rows are shared with [`crate::orbital::mean_orbital_elements`]
//! and use the Simon et al. (1994)/Meeus Table 31.A polynomials. The module
//! provides a canonical coefficient serialization, its SHA-256 identity, physical
//! vectors and J2000-anchored element recovery. No independent coefficient table
//! or ephemeris acquisition is introduced.
//!
//! Epochs are TT Julian days with `t = (jd_tt − J2000)/36525`. Distances are AU,
//! eccentricities dimensionless and public angles degrees. Source elements refer
//! to the mean ecliptic/equinox of date. At J2000 the relative reference frame is
//! the identity.
//!
//! Recovery uses binary64 geometry and is not an arbitrary-precision or interval
//! certificate. Unknown bodies, nonfinite epochs, degenerate norms and unresolved
//! angular origins return `None`. These structural reductions retain the limits
//! of the underlying secular element model.

use sha2::{Digest, Sha256};

use crate::constants;
use crate::download::hex_of;
use crate::orbital::{self, J2000_JD};

/// TT Julian centuries per Julian century (inline `36525.0` century literal in the model evaluator).
pub const CENTURY_DAYS: f64 = 36525.0;

/// Registered perturbing planets in canonical order
pub const PERTURBER_IDS: [i32; 4] = [
    constants::JUPITER,
    constants::SATURN,
    constants::URANUS,
    constants::NEPTUNE,
];

/// Canonical field order.
const CANONICAL_FIELDS: [&str; 5] = ["a", "e", "i", "node", "perihelion"];

/// Pinned SHA-256 of the canonical payload.
pub const CANONICAL_SHA256: &str =
    "f129d9b73f27d8d9bad8a81e54bd4eb25e375235039fd543e1f3d11402798bf5";

/// Canonical payload length in bytes (same pinned test).
pub const CANONICAL_LEN: usize = 852;

/// Whether `body_id` is one of the four registered perturbing planets
pub fn is_perturbing_planet(body_id: i32) -> bool {
    PERTURBER_IDS.contains(&body_id)
}

/// Registered coefficient row for one perturber (one table row: the slices are the immutable values held by the shared the model
/// table — this accessor keeps no second table). Fields are private with accessors per the validated-invariant rule;
/// every row has a non-empty polynomial per field and a present node.
pub struct PerturberRow {
    body_id: i32,
    a: &'static [f64],
    e: &'static [f64],
    i: &'static [f64],
    node: &'static [f64],
    perihelion: &'static [f64],
}

impl PerturberRow {
    /// the model table key (planet ID) for this row.
    pub fn body_id(&self) -> i32 {
        self.body_id
    }
    /// `a` field (AU polynomial, ascending powers).
    pub fn a(&self) -> &'static [f64] {
        self.a
    }
    /// `e` field (eccentricity polynomial, ascending powers).
    /// Dimensionless.
    pub fn e(&self) -> &'static [f64] {
        self.e
    }
    /// `i` field (inclination polynomial in degrees, ascending powers).
    pub fn i(&self) -> &'static [f64] {
        self.i
    }
    /// `node` field (ascending-node polynomial in degrees).
    pub fn node(&self) -> &'static [f64] {
        self.node
    }
    /// `perihelion` field (perihelion-longitude polynomial in degrees).
    pub fn perihelion(&self) -> &'static [f64] {
        self.perihelion
    }
}

/// Registered row for a perturber, or `None` for any other body
pub fn perturber_coefficients(body_id: i32) -> Option<PerturberRow> {
    if !is_perturbing_planet(body_id) {
        return None;
    }
    let el = orbital::mean_element_coeffs(body_id)?;
    let node = el.omega?;
    Some(PerturberRow {
        body_id,
        a: el.a,
        e: el.e,
        i: el.i,
        node,
        perihelion: el.pi,
    })
}

/// Canonical serialization:
/// per body, big-endian `>q` identifier, then per field the big-endian `>I`
/// name length, the UTF-8 name, the big-endian `>I` value count, and each
/// value as big-endian `>d`.
pub fn canonical_coefficients_bytes() -> Vec<u8> {
    let mut out = Vec::with_capacity(CANONICAL_LEN);
    for body_id in PERTURBER_IDS {
        let row = perturber_coefficients(body_id).expect("registered perturber row");
        out.extend_from_slice(&(body_id as i64).to_be_bytes());
        for name in CANONICAL_FIELDS {
            let values = match name {
                "a" => row.a,
                "e" => row.e,
                "i" => row.i,
                "node" => row.node,
                "perihelion" => row.perihelion,
                _ => unreachable!("canonical field table"),
            };
            out.extend_from_slice(&(name.len() as u32).to_be_bytes());
            out.extend_from_slice(name.as_bytes());
            out.extend_from_slice(&(values.len() as u32).to_be_bytes());
            for v in values {
                out.extend_from_slice(&v.to_be_bytes());
            }
        }
    }
    out
}

/// SHA-256 hex digest of the canonical payload
/// .
pub fn coefficient_source_sha256() -> String {
    hex_of(&Sha256::digest(canonical_coefficients_bytes()))
}

/// Physical source vectors at `jd_tt` (uncertified `f64` analogue of ): semi-major axis plus the date-frame eccentricity and
/// orbit-pole vectors, with the calculation order and association
/// (Horner evaluation, `* (pi / 180)` radian scaling, `(perihelion - node)`
/// argument of perihelion, `e * p` eccentricity scaling). Returns `None` for
/// non-perturbers and non-finite epochs.
pub fn source_vectors(body_id: i32, jd_tt: f64) -> Option<(f64, [f64; 3], [f64; 3])> {
    let row = perturber_coefficients(body_id)?;
    if !jd_tt.is_finite() {
        return None;
    }
    let t = (jd_tt - J2000_JD) / CENTURY_DAYS;
    let a = orbital::poly_asc(row.a, t);
    let ecc = orbital::poly_asc(row.e, t);
    let incl = orbital::poly_asc(row.i, t) * (std::f64::consts::PI / 180.0);
    let node = orbital::poly_asc(row.node, t) * (std::f64::consts::PI / 180.0);
    let peri = (orbital::poly_asc(row.perihelion, t) - orbital::poly_asc(row.node, t))
        * (std::f64::consts::PI / 180.0);
    let (sin_node, cos_node) = node.sin_cos();
    let (sin_i, cos_i) = incl.sin_cos();
    let (sin_peri, cos_peri) = peri.sin_cos();
    let px = cos_node * cos_peri - sin_node * cos_i * sin_peri;
    let py = sin_node * cos_peri + cos_node * cos_i * sin_peri;
    let pz = sin_i * sin_peri;
    let pole = [sin_i * sin_node, -sin_i * cos_node, cos_i];
    let evec = [ecc * px, ecc * py, ecc * pz];
    Some((a, evec, pole))
}

/// Degeneracy proof class of a recovered tuple. For the registered perturbers at finite dates the class is
/// always [`ProofClass::Ordinary`]; the canonical branches are reachable only
/// through constructed degenerate inputs, and the f64 zero predicates below
/// are representational checks, not exact-ball certificates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProofClass {
    /// Exact retrograde pole (`-1`) with exact-zero eccentricity-plane radius.
    CanonicalPole,
    /// Exact equatorial pole with exact-zero eccentricity-plane radius.
    CanonicalEquator,
    /// Exact-zero eccentricity vector on an inclined orbit.
    ZeroEccentricity,
    /// The generic inclined, eccentric case.
    Ordinary,
}

impl ProofClass {
    /// `proof_class` string.
    pub fn as_str(self) -> &'static str {
        match self {
            ProofClass::CanonicalPole => "canonical-pole",
            ProofClass::CanonicalEquator => "canonical-equator",
            ProofClass::ZeroEccentricity => "zero-eccentricity",
            ProofClass::Ordinary => "ordinary",
        }
    }
}

/// Euclidean norm with the grouping
/// (`(x*x + y*y) + z*z`).
fn norm3(v: [f64; 3]) -> f64 {
    ((v[0] * v[0] + v[1] * v[1]) + v[2] * v[2]).sqrt()
}

/// Normalize the heliocentric mean longitude with Euclidean modulo.
///
/// This binary64 geometry is used for structural bounds and diagnostics;
/// it is not an interval-arithmetic proof or numerical enclosure.
fn angle_degrees(pair_x: f64, pair_y: f64) -> Option<f64> {
    if pair_x == 0.0 {
        if pair_y == 0.0 {
            return None;
        }
        return Some(if pair_y > 0.0 { 90.0 } else { 270.0 });
    }
    if pair_y == 0.0 {
        return Some(if pair_x > 0.0 { 0.0 } else { 180.0 });
    }
    let deg = pair_y.atan2(pair_x).to_degrees().rem_euclid(360.0);
    deg.is_finite().then_some(deg)
}

/// Recover `(a_au, e, i_deg, node_deg, omega_deg)` from J2000-anchored vectors
/// (uncertified `f64` analogue of ): same branch
/// structure — pole normalization with a positive-finite gate, `acos`-domain
/// gate, canonical-equator/pole fast paths, node-basis/quadrature projection
/// for the generic case. Returns `None` where the model raises
/// `CalculationError` (degenerate norms, non-finite eccentricity, pole outside
/// the `acos` domain, unresolved angle origin).
pub fn recover_elements(
    a: f64,
    evec: [f64; 3],
    pole_value: [f64; 3],
) -> Option<([f64; 5], ProofClass)> {
    if !a.is_finite() {
        return None;
    }
    let pole_norm = norm3(pole_value);
    if !pole_norm.is_finite() || pole_norm <= 0.0 {
        return None;
    }
    let pole = [
        pole_value[0] / pole_norm,
        pole_value[1] / pole_norm,
        pole_value[2] / pole_norm,
    ];
    let ecc = norm3(evec);
    if !ecc.is_finite() || ecc < 0.0 {
        return None;
    }
    let pole_z = pole[2];
    if pole_z < -1.0 || pole_z > 1.0 {
        return None;
    }
    let inclination = pole_z.acos();
    let sin_i = (pole[0] * pole[0] + pole[1] * pole[1]).sqrt();
    // Representational zeros only: the model proves exact-ball zeros here.
    let exact_zero_e = ecc == 0.0;
    let exact_zero_i = sin_i == 0.0;
    let exact_retrograde = exact_zero_i && pole[2] == -1.0;
    if exact_zero_i {
        let omega_deg = if exact_zero_e {
            0.0
        } else if exact_retrograde {
            angle_degrees(evec[0], -evec[1])?
        } else {
            angle_degrees(evec[0], evec[1])?
        };
        let tuple = (
            a,
            ecc,
            if exact_retrograde { 180.0 } else { 0.0 },
            0.0,
            omega_deg,
        );
        let class = if exact_retrograde {
            ProofClass::CanonicalPole
        } else {
            ProofClass::CanonicalEquator
        };
        return Some(([tuple.0, tuple.1, tuple.2, tuple.3, tuple.4], class));
    }
    let node_deg = angle_degrees(-pole[1], pole[0])?;
    let node_basis = [-pole[1] / sin_i, pole[0] / sin_i, 0.0];
    let quadrature = [
        pole[1] * node_basis[2] - pole[2] * node_basis[1],
        pole[2] * node_basis[0] - pole[0] * node_basis[2],
        pole[0] * node_basis[1] - pole[1] * node_basis[0],
    ];
    let omega_deg = if exact_zero_e {
        0.0
    } else {
        let peri = [evec[0] / ecc, evec[1] / ecc, evec[2] / ecc];
        let dot_basis =
            (peri[0] * node_basis[0] + peri[1] * node_basis[1]) + peri[2] * node_basis[2];
        let dot_quad =
            (peri[0] * quadrature[0] + peri[1] * quadrature[1]) + peri[2] * quadrature[2];
        angle_degrees(dot_basis, dot_quad)?
    };
    let tuple = [
        a,
        ecc,
        inclination * (180.0 / std::f64::consts::PI),
        node_deg,
        omega_deg,
    ];
    let class = if exact_zero_e {
        ProofClass::ZeroEccentricity
    } else {
        ProofClass::Ordinary
    };
    Some((tuple, class))
}

/// Anchored-J2000 tuple for one perturber (uncertified `f64` analogue of the
/// J2000 path): at `jd_tt == J2000` the
/// returns the literal identity, so the date-frame
/// source vectors are already J2000-anchored and only `_recover_tuple`
/// remains. Returns `None` for non-perturbers and non-finite epochs; any other
/// date needs the ERFA relative frame and a separate relative-frame construction; this helper is J2000-only.
pub fn anchored_tuple_j2000(body_id: i32) -> Option<([f64; 5], ProofClass)> {
    let (a, evec, pole) = source_vectors(body_id, J2000_JD)?;
    recover_elements(a, evec, pole)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// B-11 registered rows, transcribed from the pinned V 4 focus test
    /// (`tests/test_perturber_elements.py` @ `e6d7b18...`); the values are
    /// value-identical to the rows, which is itself
    /// asserted by the agreement of `perturber_coefficients` with these
    /// literals below.
    const EXPECTED: [(i32, [&[f64]; 5]); 4] = [
        (
            constants::JUPITER,
            [
                &[5.202603191, 0.0000001913],
                &[0.04849485, 0.000163244, -0.0000004719, -0.00000000197],
                &[1.303270, -0.0054966, 0.00000465, -0.000000004],
                &[100.464441, 1.0209550, 0.00040117, 0.000000569],
                &[14.331309, 1.6126668, 0.00103127, -0.000004569],
            ],
        ),
        (
            constants::SATURN,
            [
                &[9.554909596, -0.0000021389],
                &[0.05550862, -0.000346818, -0.0000006456, 0.00000000338],
                &[2.488878, -0.0037363, -0.00001516, 0.000000089],
                &[113.665524, 0.8770979, -0.00012067, -0.000002380],
                &[93.056787, 1.9637694, 0.00083757, 0.000004899],
            ],
        ),
        (
            constants::URANUS,
            [
                &[19.218446062, -0.0000000372, 0.00000000098],
                &[0.04629590, -0.000027337, 0.0000000790, 0.00000000025],
                &[0.773196, 0.0007744, 0.00003749, -0.000000092],
                &[74.005947, 0.5211258, 0.00133982, 0.000018516],
                &[173.005159, 1.4863784, 0.00021450, 0.000000433],
            ],
        ),
        (
            constants::NEPTUNE,
            [
                &[30.110386869, -0.0000001663, 0.00000000069],
                &[0.00898809, 0.000006408, -0.0000000008, -0.00000000005],
                &[1.769952, -0.0093082, -0.00000708, 0.000000028],
                &[131.784057, 1.1022057, 0.00026006, -0.000000636],
                &[48.123691, 1.4262677, 0.00037918, -0.000000003],
            ],
        ),
    ];

    #[test]
    fn accessor_is_source_identical() {
        for (body_id, fields) in EXPECTED {
            assert!(is_perturbing_planet(body_id));
            let row = perturber_coefficients(body_id).expect("registered perturber");
            assert_eq!(row.body_id(), body_id);
            assert_eq!(row.a(), fields[0]);
            assert_eq!(row.e(), fields[1]);
            assert_eq!(row.i(), fields[2]);
            assert_eq!(row.node(), fields[3]);
            assert_eq!(row.perihelion(), fields[4]);
        }
        assert_eq!(PERTURBER_IDS.len(), 4);
    }

    #[test]
    fn non_perturbers_and_bad_epochs_are_none() {
        for body_id in [
            constants::MERCURY,
            constants::EARTH,
            constants::MARS,
            constants::PLUTO,
            0,
            -1,
        ] {
            assert!(!is_perturbing_planet(body_id));
            assert!(perturber_coefficients(body_id).is_none());
            assert!(source_vectors(body_id, J2000_JD).is_none());
            assert!(anchored_tuple_j2000(body_id).is_none());
        }
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(source_vectors(constants::JUPITER, bad).is_none());
        }
    }

    #[test]
    fn canonical_encoding_shape_and_order() {
        let payload = canonical_coefficients_bytes();
        assert_eq!(payload.len(), CANONICAL_LEN);
        // Body order is Jupiter, Saturn, Uranus, Neptune.
        assert_eq!(&payload[0..8], &(constants::JUPITER as i64).to_be_bytes());
        // Fully decode: every block must round-trip its registered row.
        let mut offset = 0usize;
        for body_id in PERTURBER_IDS {
            let id = i64::from_be_bytes(payload[offset..offset + 8].try_into().unwrap());
            assert_eq!(id as i32, body_id);
            offset += 8;
            let row = perturber_coefficients(body_id).expect("registered perturber");
            let fields = [
                ("a", row.a()),
                ("e", row.e()),
                ("i", row.i()),
                ("node", row.node()),
                ("perihelion", row.perihelion()),
            ];
            for (name, values) in fields {
                let name_len =
                    u32::from_be_bytes(payload[offset..offset + 4].try_into().unwrap()) as usize;
                offset += 4;
                assert_eq!(&payload[offset..offset + name_len], name.as_bytes());
                offset += name_len;
                let count =
                    u32::from_be_bytes(payload[offset..offset + 4].try_into().unwrap()) as usize;
                offset += 4;
                assert_eq!(count, values.len());
                for expected in values {
                    let word = f64::from_be_bytes(payload[offset..offset + 8].try_into().unwrap());
                    assert_eq!(word.to_bits(), expected.to_bits());
                    offset += 8;
                }
            }
        }
        assert_eq!(offset, CANONICAL_LEN);
    }

    #[test]
    fn canonical_digest_pin() {
        assert_eq!(coefficient_source_sha256(), CANONICAL_SHA256);
    }

    /// A-priori agreement bounds for the uncertified J2000 analogue vs the
    /// pinned B-11 certificate `(5.202603191, 0.04849485, 1.30327, 100.464441,
    /// 273.866868)`: `a` is an exact Horner constant (exactly); `e` crosses
    /// one sqrt of a near-unit vector (~2 ulps, bound 1 e-15 absolute); the
    /// three angles each cross one sin/cos pair, one `atan2`, and two radian
    /// scalings (~10 ulps worst case, bound 1 e-12 deg absolute). These bounds
    /// are declared from the operation count before running, not fitted.
    #[test]
    fn j2000_jupiter_matches_certificate() {
        let (tuple, class) = anchored_tuple_j2000(constants::JUPITER).expect("Jupiter row");
        assert_eq!(class, ProofClass::Ordinary);
        assert_eq!(tuple[0].to_bits(), 5.202603191f64.to_bits());
        assert!(
            (tuple[1] - 0.04849485).abs() <= 1e-15,
            "e slot: {}",
            tuple[1]
        );
        assert!((tuple[2] - 1.30327).abs() <= 1e-12, "i slot: {}", tuple[2]);
        assert!(
            (tuple[3] - 100.464441).abs() <= 1e-12,
            "node slot: {}",
            tuple[3]
        );
        assert!(
            (tuple[4] - 273.866868).abs() <= 1e-12,
            "omega slot: {}",
            tuple[4]
        );
    }

    #[test]
    fn j2000_all_perturbers_ordinary_and_finite() {
        for body_id in PERTURBER_IDS {
            let (tuple, class) = anchored_tuple_j2000(body_id).expect("perturber row");
            assert_eq!(class, ProofClass::Ordinary, "body {body_id}");
            assert!(tuple.iter().all(|v| v.is_finite()), "body {body_id}");
            // Semi-major axis is the exact constant term.
            let row = perturber_coefficients(body_id).expect("perturber row");
            assert_eq!(tuple[0].to_bits(), row.a()[0].to_bits());
        }
    }

    #[test]
    fn j2000_source_vectors_match_constants() {
        let (a, evec, pole) = source_vectors(constants::JUPITER, J2000_JD).expect("row");
        assert_eq!(a.to_bits(), 5.202603191f64.to_bits());
        // The date-frame vectors are unit-scaled by construction: |pole| is
        // sin^2+cos^2-grade (~2 ulps) and |evec| is e times a unit p-vector.
        assert!((norm3(pole) - 1.0).abs() <= 1e-15);
        assert!((norm3(evec) - 0.04849485).abs() <= 1e-15);
    }

    #[test]
    fn recover_refuses_degenerate_inputs() {
        // Non-finite semi-major axis, zero pole vector, pole norm overflowing
        // to infinity, and a non-finite pole component. (A finite nonzero
        // pole always normalizes inside the acos domain, so the domain gate
        // mirrors the model structure without a reachable finite trigger.)
        assert!(recover_elements(1.0, [0.01, 0.0, 0.0], [0.0, 0.0, 0.0]).is_none());
        assert!(recover_elements(1.0, [0.01, 0.0, 0.0], [1.0, 0.0, 1e300]).is_none());
        assert!(recover_elements(f64::NAN, [0.01, 0.0, 0.0], [0.0, 0.0, 1.0]).is_none());
        assert!(recover_elements(1.0, [0.01, 0.0, 0.0], [f64::NAN, 0.0, 1.0]).is_none());
    }

    #[test]
    fn recover_canonical_branches() {
        // Exact equatorial pole, eccentric: omega from the evec longitude.
        // (`e` crosses one sqrt, so it is pinned within 3 ulps, not exactly.)
        let (tuple, class) =
            recover_elements(2.5, [0.03, 0.0, 0.0], [0.0, 0.0, 1.0]).expect("canonical");
        assert_eq!(class, ProofClass::CanonicalEquator);
        assert_eq!(tuple[0].to_bits(), 2.5f64.to_bits());
        assert!((tuple[1] - 0.03).abs() <= 1e-17, "e slot: {}", tuple[1]);
        assert_eq!([tuple[2], tuple[3], tuple[4]], [0.0, 0.0, 0.0]);
        // Exact retrograde pole, eccentric: omega from the mirrored longitude.
        let (tuple, class) =
            recover_elements(2.5, [0.0, 0.03, 0.0], [0.0, 0.0, -1.0]).expect("retrograde");
        assert_eq!(class, ProofClass::CanonicalPole);
        assert_eq!(tuple[2], 180.0);
        assert!((tuple[4] - 270.0).abs() <= 1e-12, "omega: {}", tuple[4]);
        // Inclined orbit, exact-zero eccentricity: zero-eccentricity class.
        let pole = [0.0, -0.02, 0.99979998];
        let (tuple, class) = recover_elements(1.0, [0.0, 0.0, 0.0], pole).expect("zero-e");
        assert_eq!(class, ProofClass::ZeroEccentricity);
        assert_eq!(tuple[4], 0.0);
    }

    #[test]
    fn proof_class_spellings() {
        assert_eq!(ProofClass::CanonicalPole.as_str(), "canonical-pole");
        assert_eq!(ProofClass::CanonicalEquator.as_str(), "canonical-equator");
        assert_eq!(ProofClass::ZeroEccentricity.as_str(), "zero-eccentricity");
        assert_eq!(ProofClass::Ordinary.as_str(), "ordinary");
    }
}
