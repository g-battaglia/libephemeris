//! Curated minor-body elements, Keplerian propagation and discovery helpers.
//!
//! The model combines curated and multi-epoch osculating elements, epoch blending,
//! secular perturbation rates and native Vondrak frame transformations. Kernel
//! positions, when supplied, precede the approximate Keplerian model according
//! to the caller's serving policy. Network acquisition remains an explicit,
//! policy-controlled operation through the native transport.
//!
//! Epochs are TT Julian days, distances AU and angles degrees; Cartesian model
//! positions use heliocentric J2000 ecliptic axes unless stated otherwise.
//! Nutation in longitude is an explicit input to the of-date composition.
//! Approximate main-belt accuracy is tens of arcseconds and distant-object
//! accuracy is arcminutes near the element epochs; errors grow off epoch and
//! resonant phases need their own model.
//!
//! Unknown elements and unavailable optional discoveries return `None` where
//! documented. Nonfinite inputs, invalid orbital domains, strict provider
//! failures and sealed-network refusals retain typed errors. These models do
//! not silently supply the precision of an unavailable numerical integration.

use crate::constants::{
    APOPHIS, BENNU, CERES, CHIRON, EROS, GONGGONG, ITOKAWA, IXION, JUNO, ORCUS, PALLAS, RYUGU,
    TOUTATIS, VESTA,
};
use crate::errors::Error;

type Result<T> = std::result::Result<T, Error>;

/// Provider queried before the Keplerian fallback: returns the served
/// heliocentric `(lon_deg, lat_deg, dist_au)` triple, `Ok(None)` when the
/// body is not registered, or a range-type miss that likewise falls
/// through (every other failure propagates).
pub type SpkPositionProvider<'a> = &'a dyn Fn(i32, f64) -> Result<Option<(f64, f64, f64)>>;

/// SBDB element lookup by catalog number: `Ok(None)` for unknown bodies.
pub type SbdbElementsProvider<'a> = &'a dyn Fn(i64) -> Result<Option<OrbitalElements>>;

/// SBDB name-to-number lookup over the raw (unnormalized) name.
pub type AsteroidNameProvider<'a> = &'a dyn Fn(&str) -> Result<Option<i64>>;

// ---------------------------------------------------------------------------
// Physical constants and planet data
// ---------------------------------------------------------------------------

/// Gaussian `k^2` solar mass parameter, AU^3/day^2.
pub const GM_SUN: f64 = 0.00029591220828559;
/// Jupiter/Sun mass ratio.
pub const MASS_RATIO_JUPITER: f64 = 1.0 / 1047.348644;
/// Saturn/Sun mass ratio.
pub const MASS_RATIO_SATURN: f64 = 1.0 / 3497.901768;
/// Uranus/Sun mass ratio.
pub const MASS_RATIO_URANUS: f64 = 1.0 / 22902.0;
/// Neptune/Sun mass ratio. Mass ratios follow the JPL DE440/DE441 planetary
/// masses (Park, Folkner, Williams & Boggs 2021).
pub const MASS_RATIO_NEPTUNE: f64 = 1.0 / 19412.0;

/// J2000.0 mean elements + Simon et al. (1994) linear rates for the four
/// giant planets: `(a_AU, mass_ratio, e, omega_deg, node_deg, i_deg)`.
/// `JUPITER_OMEGA`/`SATURN_OMEGA` carry the model "mean argument of
/// perihelion" vocabulary (varpi = omega + node in the comments there).
/// Mean longitudes of perihelion are 14.7 (Jupiter), 93.1 (Saturn), 170.9
/// (Uranus) and 44.9 (Neptune) degrees; the element sources follow the
/// Standish (1992/2006) mean-element fits.
pub const JUPITER_A: f64 = 5.2026;
/// Jupiter J2000.0 mean eccentricity.
pub const JUPITER_E: f64 = 0.0485;
/// Jupiter J2000.0 mean inclination, degrees.
pub const JUPITER_I: f64 = 1.303;
/// Jupiter mean argument of perihelion, degrees.
pub const JUPITER_OMEGA: f64 = 274.25;
/// Jupiter mean longitude of ascending node, degrees.
pub const JUPITER_NODE: f64 = 100.46;
/// Jupiter mean motion, degrees/day.
pub const JUPITER_N: f64 = 0.08309;
/// Saturn J2000.0 mean semi-major axis, AU (Standish fit).
pub const SATURN_A: f64 = 9.5549;
/// Saturn J2000.0 mean eccentricity.
pub const SATURN_E: f64 = 0.0555;
/// Saturn J2000.0 mean inclination, degrees.
pub const SATURN_I: f64 = 2.489;
/// Saturn mean argument of perihelion, degrees.
pub const SATURN_OMEGA: f64 = 339.39;
/// Saturn mean longitude of ascending node, degrees.
pub const SATURN_NODE: f64 = 113.66;
/// Saturn mean motion, degrees/day.
pub const SATURN_N: f64 = 0.03346;
/// Uranus J2000.0 mean semi-major axis, AU (Standish fit).
pub const URANUS_A: f64 = 19.2184;
/// Uranus J2000.0 mean eccentricity.
pub const URANUS_E: f64 = 0.0457;
/// Uranus J2000.0 mean inclination, degrees.
pub const URANUS_I: f64 = 0.772;
/// Uranus mean argument of perihelion, degrees.
pub const URANUS_OMEGA: f64 = 96.9;
/// Uranus mean longitude of ascending node, degrees.
pub const URANUS_NODE: f64 = 74.0;
/// Uranus mean motion, degrees/day.
pub const URANUS_N: f64 = 0.01177;
/// Neptune J2000.0 mean semi-major axis, AU (Standish fit).
pub const NEPTUNE_A: f64 = 30.1104;
/// Neptune J2000.0 mean eccentricity.
pub const NEPTUNE_E: f64 = 0.0086;
/// Neptune J2000.0 mean inclination, degrees.
pub const NEPTUNE_I: f64 = 1.769;
/// Neptune mean argument of perihelion, degrees.
pub const NEPTUNE_OMEGA: f64 = 273.2;
/// Neptune mean longitude of ascending node, degrees.
pub const NEPTUNE_NODE: f64 = 131.7;
/// Neptune mean motion, degrees/day.
pub const NEPTUNE_N: f64 = 0.006021;

/// Simon et al. (1994) A&A 282, 663 century rates with Standish mean
/// elements, paired with the elements above:
/// `(de, di_deg, domega_deg, dnode_deg)` per Julian century (36525 days).
pub const JUPITER_DE: f64 = -0.00036;
/// Jupiter inclination rate, degrees per Julian century.
pub const JUPITER_DI: f64 = -0.0019;
/// Jupiter longitude-of-perihelion rate, degrees per Julian century.
pub const JUPITER_DOMEGA: f64 = 0.3625;
/// Jupiter ascending-node rate, degrees per Julian century.
pub const JUPITER_DNODE: f64 = 0.3574;
/// Saturn eccentricity rate per Julian century (dimensionless/century).
pub const SATURN_DE: f64 = -0.00051;
/// Saturn inclination rate, degrees per Julian century.
pub const SATURN_DI: f64 = 0.0019;
/// Saturn longitude-of-perihelion rate, degrees per Julian century.
pub const SATURN_DOMEGA: f64 = 0.7375;
/// Saturn ascending-node rate, degrees per Julian century.
pub const SATURN_DNODE: f64 = -0.2510;
/// Uranus eccentricity rate per Julian century (dimensionless/century).
pub const URANUS_DE: f64 = -0.00018;
/// Uranus inclination rate, degrees per Julian century.
pub const URANUS_DI: f64 = -0.0024;
/// Uranus longitude-of-perihelion rate, degrees per Julian century.
pub const URANUS_DOMEGA: f64 = 0.0710;
/// Uranus ascending-node rate, degrees per Julian century.
pub const URANUS_DNODE: f64 = 0.0174;
/// Neptune eccentricity rate per Julian century (dimensionless/century).
pub const NEPTUNE_DE: f64 = 0.00007;
/// Neptune inclination rate, degrees per Julian century.
pub const NEPTUNE_DI: f64 = -0.0035;
/// Neptune longitude-of-perihelion rate, degrees per Julian century.
pub const NEPTUNE_DOMEGA: f64 = -0.0070;
/// Neptune ascending-node rate, degrees per Julian century.
pub const NEPTUNE_DNODE: f64 = -0.0063;

/// J2000.0 reference epoch for the libration/Neptune models.
pub const J2000_EPOCH: f64 = 2451545.0;
/// Neptune mean longitude at J2000.0, degrees (Standish fit).
pub const NEPTUNE_MEAN_LONGITUDE_J2000: f64 = 304.88003;
/// Neptune mean motion, degrees/day (resonance reference).
pub const NEPTUNE_MEAN_MOTION: f64 = 0.006021;
/// Default fractional tolerance for resonance detection.
pub const DEFAULT_RESONANCE_TOLERANCE: f64 = 0.02;

/// Resonant semi-major axis for a `p:q` resonance with Neptune: `a_neptune * (q / p)^(2/3)`.
pub fn resonant_semi_major_axis(p: i32, q: i32) -> f64 {
    NEPTUNE_A * ((q as f64) / (p as f64)).powf(2.0 / 3.0)
}

/// Known Neptune mean-motion resonances: `(p, q, name, tolerance)`.
/// The resonant axis follows from [`resonant_semi_major_axis`].
pub const NEPTUNE_RESONANCES: [(i32, i32, &str, f64); 7] = [
    (3, 4, "3:4 resonance", 0.02),
    (3, 5, "3:5 resonance", 0.02),
    (2, 3, "plutino", 0.02),
    (4, 7, "4:7 resonance", 0.02),
    (1, 2, "twotino", 0.02),
    (2, 5, "2:5 resonance", 0.02),
    (1, 3, "1:3 resonance", 0.02),
];

/// Libration parameters for a resonant TNO.
#[derive(Clone, Copy, Debug)]
pub struct LibrationParameters {
    /// Maximum angular displacement from center, degrees.
    pub amplitude: f64,
    /// Libration period, days.
    pub period: f64,
    /// Center value of the resonant argument, degrees.
    pub center: f64,
    /// Phase at J2000.0, radians.
    pub phase_j2000: f64,
    /// `p` in the `p:q` resonance.
    pub resonance_p: i32,
    /// `q` in the `p:q` resonance.
    pub resonance_q: i32,
}

/// Project libration fits for known resonant TNOs. Gonggong's row is a
/// documented project estimate, not a published fit.
pub fn libration_parameters(body_id: i32) -> Option<LibrationParameters> {
    match body_id {
        x if x == IXION => Some(LibrationParameters {
            amplitude: 78.0,
            period: 7_233_000.0,
            center: 180.0,
            phase_j2000: 2.14,
            resonance_p: 2,
            resonance_q: 3,
        }),
        x if x == ORCUS => Some(LibrationParameters {
            amplitude: 68.0,
            period: 7_379_000.0,
            center: 180.0,
            phase_j2000: 4.71,
            resonance_p: 2,
            resonance_q: 3,
        }),
        x if x == GONGGONG => Some(LibrationParameters {
            amplitude: 30.0,
            period: 9_131_000.0,
            center: 180.0,
            phase_j2000: 1.05,
            resonance_p: 3,
            resonance_q: 10,
        }),
        _ => None,
    }
}

/// Whether a body carries resonant libration parameters.
pub fn has_libration_model(body_id: i32) -> bool {
    libration_parameters(body_id).is_some()
}

/// Libration parameters for a resonant body, if any.
pub fn get_libration_parameters(body_id: i32) -> Option<LibrationParameters> {
    libration_parameters(body_id)
}

/// Measured Yarkovsky `da/dt` in AU per million years. Per-row basis: Apophis
/// (Brozovic et al. 2018), Bennu (OSIRIS-REx), Ryugu (Hayabusa 2 estimate),
/// Itokawa (Vokrouhlicky et al. 2008); Eros and Toutatis are estimates.
pub fn yarkovsky_da_dt_au_per_my(body_id: i32) -> Option<f64> {
    match body_id {
        x if x == APOPHIS => Some(-2.899e-4),
        x if x == BENNU => Some(-18.99e-4),
        x if x == RYUGU => Some(-3.5e-4),
        x if x == ITOKAWA => Some(-3.5e-4),
        x if x == EROS => Some(-0.5e-4),
        x if x == TOUTATIS => Some(-1.5e-4),
        _ => None,
    }
}

// ---------------------------------------------------------------------------
// Orbital elements
// ---------------------------------------------------------------------------

/// Classical Keplerian orbital elements for a minor body: osculating elements at [`OrbitalElements::epoch`].
///
/// Angles in degrees, `a` in AU (`q`, the perihelion distance, for
/// parabolic rows), `n` in degrees/day.
#[derive(Clone, Debug)]
pub struct OrbitalElements {
    name: String,
    epoch: f64,
    a: f64,
    e: f64,
    i: f64,
    omega: f64,
    node: f64,
    m0: f64,
    n: f64,
    // Set only by an actual catalogue/body acquisition, never by the public
    // scalar constructors. This identity keeps a selected asteroid's own
    // force out of its IVP without mistaking a close encounter for self-force.
    physical_spice_id: Option<i32>,
}

impl OrbitalElements {
    /// Checked constructor: every numeric field must be finite.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        name: String,
        epoch: f64,
        a: f64,
        e: f64,
        i: f64,
        omega: f64,
        node: f64,
        m0: f64,
        n: f64,
    ) -> Result<Self> {
        for v in [epoch, a, e, i, omega, node, m0, n] {
            if !v.is_finite() {
                return Err(Error::NonFiniteArgument { value: v });
            }
        }
        Ok(Self {
            name,
            epoch,
            a,
            e,
            i,
            omega,
            node,
            m0,
            n,
            physical_spice_id: None,
        })
    }

    /// Own the source N-body dataclass fields without numerical admission.
    ///
    /// Epoch is scalar TT JD; a is AU; e is dimensionless; i/omega/node/m 0 are
    /// degrees; n is the original stored degrees/day field, never derived here.
    /// Every IEEE value, including signed zero, NaN and infinity, is retained.
    /// Names are owned UTF-8, including embedded NUL. Native timestep/hash/orbit
    /// initialization owns the corresponding source conversion errors; no
    /// validation or physical-sanity advisory is executed by construction.
    ///
    /// This explicit interoperability record is separate from [`Self::new`]'s
    /// unchanged finite-field contract. Other kernels keep their own checked
    /// domains. Constructing a record does not certify a physically valid orbit
    /// or that any native integration will finish. Caller retains no borrow.
    #[allow(clippy::too_many_arguments)]
    pub fn from_nbody_source_record(
        name: String,
        epoch: f64,
        a: f64,
        e: f64,
        i: f64,
        omega: f64,
        node: f64,
        m0: f64,
        n: f64,
    ) -> Self {
        Self {
            name,
            epoch,
            a,
            e,
            i,
            omega,
            node,
            m0,
            n,
            physical_spice_id: None,
        }
    }

    /// Stored body name (see [`OrbitalElements`] for the element contract).
    pub fn name(&self) -> &str {
        &self.name
    }
    /// Osculating epoch as a TT Julian day; anomalies propagate from here.
    pub fn epoch(&self) -> f64 {
        self.epoch
    }
    /// Semi-major axis in AU (perihelion distance `q` for parabolic rows).
    pub fn a(&self) -> f64 {
        self.a
    }
    /// Eccentricity (dimensionless).
    pub fn e(&self) -> f64 {
        self.e
    }
    /// Inclination in degrees.
    pub fn i(&self) -> f64 {
        self.i
    }
    /// Argument of perihelion in degrees (`varpi = omega + node`).
    pub fn omega(&self) -> f64 {
        self.omega
    }
    /// Longitude of ascending node.
    pub fn node(&self) -> f64 {
        self.node
    }
    /// Mean anomaly at the osculating epoch in degrees.
    pub fn m0(&self) -> f64 {
        self.m0
    }
    /// Mean daily motion in degrees per day.
    pub fn n(&self) -> f64 {
        self.n
    }

    /// Physical target retained by catalogue acquisition, without inference
    /// from the record's name, numerical fields or proximity to a perturber.
    pub(crate) fn physical_spice_id(&self) -> Option<i32> {
        self.physical_spice_id
    }

    /// Restore only the identity already retained by the internal Python
    /// bridge for this owned catalogue object; scalar constructors stay anonymous.
    #[cfg(feature = "python-internals")]
    pub(crate) fn with_retained_physical_identity(mut self, target: Option<i32>) -> Self {
        self.physical_spice_id = target;
        self
    }

    /// Physical-sanity warnings:
    /// returns warning strings, never raises. Shapes mirror the model
    /// f-strings.
    pub fn validate(&self) -> Vec<String> {
        let mut warnings = Vec::new();
        if self.a <= 0.0 {
            warnings.push(format!(
                "{}: semi-major axis a={:.6} AU is non-positive",
                self.name, self.a
            ));
        }
        if self.e < 0.0 {
            warnings.push(format!(
                "{}: eccentricity e={:.6} is negative",
                self.name, self.e
            ));
        } else if self.e >= 1.0 && self.a > 0.0 {
            warnings.push(format!(
                "{}: eccentricity e={:.6} >= 1 but a > 0 (hyperbolic orbits should have a < 0)",
                self.name, self.e
            ));
        }
        if self.i < 0.0 || self.i > 180.0 {
            warnings.push(format!(
                "{}: inclination i={:.4}° outside [0°, 180°]",
                self.name, self.i
            ));
        }
        if self.e < 1.0 && self.n <= 0.0 {
            warnings.push(format!(
                "{}: mean motion n={:.8}°/day is non-positive for elliptic orbit",
                self.name, self.n
            ));
        }
        warnings
    }
}

/// Convenience wrapper around [`OrbitalElements::validate`].
pub fn validate_elements(elements: &OrbitalElements) -> Vec<String> {
    elements.validate()
}

/// Log physical element-validation diagnostics without changing the elements.
/// Orbital units and validation criteria are those of [`OrbitalElements`].
/// This diagnostic returns no calculated position and introduces no fallback.
pub fn validate_elements_and_log(elements: &OrbitalElements) {
    for issue in elements.validate() {
        log::warn!("Orbital element validation: {issue}");
    }
}

// ---------------------------------------------------------------------------
// Resonance detection
// ---------------------------------------------------------------------------

/// A Neptune mean-motion resonance match.
#[derive(Clone, Debug)]
pub struct ResonanceInfo {
    /// Body's period multiplier.
    pub p: i32,
    /// Neptune's period multiplier.
    pub q: i32,
    /// Common name (e.g. `"plutino"`).
    pub name: &'static str,
    /// Expected semi-major axis, AU.
    pub a_resonant: f64,
    /// Fractional detection tolerance.
    pub tolerance: f64,
}

/// Resonance detection outcome.
#[derive(Clone, Debug)]
pub struct ResonanceResult {
    /// Whether a listed Neptune resonance matched within the requested tolerance.
    pub is_resonant: bool,
    /// Closest matching catalogue resonance, or `None` when no match was found.
    pub resonance: Option<ResonanceInfo>,
    /// Fractional mean-motion-ratio deviation from the selected resonance.
    /// Zero when no resonance was selected; it is not an accuracy estimate.
    pub deviation: f64,
    /// Advisory about secular-theory limitations, present only for a match.
    pub warning_message: Option<String>,
}

/// Detect a near mean-motion resonance with Neptune: compares the measured mean-motion
/// ratio against each known `p:q` resonance within `tolerance`. Bodies
/// with `a < 20` AU are never resonant. Method: Murray & Dermott ch. 8;
/// Malhotra (1995). Secular theory degrades inside resonance (the warning
/// is returned, not raised).
pub fn detect_mean_motion_resonance(elements: &OrbitalElements, tolerance: f64) -> ResonanceResult {
    if elements.a < 20.0 {
        return ResonanceResult {
            is_resonant: false,
            resonance: None,
            deviation: 0.0,
            warning_message: None,
        };
    }
    let n_ratio = elements.n / NEPTUNE_MEAN_MOTION;
    let ratios: Vec<f64> = NEPTUNE_RESONANCES
        .iter()
        .map(|(p, q, _, _)| (*p as f64) / (*q as f64))
        .collect();
    let best =
        resonance_match_from_ratios(n_ratio, tolerance, &ratios).map(|(index, deviation)| {
            let (p, q, name, _) = NEPTUNE_RESONANCES[index];
            (
                ResonanceInfo {
                    p,
                    q,
                    name,
                    a_resonant: resonant_semi_major_axis(p, q),
                    tolerance,
                },
                deviation,
            )
        });
    match best {
        Some((resonance, deviation)) => {
            let warning = format!(
                "Body '{}' is near the {}:{} mean motion resonance with Neptune ({}). Secular perturbation theory may be inaccurate for resonant bodies. Deviation from exact resonance: {:.2}%",
                elements.name,
                resonance.p,
                resonance.q,
                resonance.name,
                deviation * 100.0
            );
            ResonanceResult {
                is_resonant: true,
                resonance: Some(resonance),
                deviation,
                warning_message: Some(warning),
            }
        }
        None => ResonanceResult {
            is_resonant: false,
            resonance: None,
            deviation: 0.0,
            warning_message: None,
        },
    }
}

/// Select a resonance after the caller has consumed the outer-axis guard.
/// `mean_motion_ratio` and each expected `p/q` are dimensionless acquired
/// values. The earliest strictly better deviation below `tolerance` wins;
/// NaNs fail the comparisons. Returns the original row index and fractional
/// deviation without constructing, copying or changing the caller's record.
pub fn resonance_match_from_ratios(
    mean_motion_ratio: f64,
    tolerance: f64,
    expected_ratios: &[f64],
) -> Option<(usize, f64)> {
    let mut best = None;
    let mut best_deviation = f64::INFINITY;
    for (index, expected) in expected_ratios.iter().copied().enumerate() {
        let deviation = (mean_motion_ratio - expected).abs() / expected;
        if deviation < tolerance && deviation < best_deviation {
            best = Some((index, deviation));
            best_deviation = deviation;
        }
    }
    best
}

/// Quick resonance check by body id:
/// unknown ids are [`Error::UnknownBody`]. The caller passes `tolerance`
/// explicitly (no default); the model default `0.02` lives in
/// `DEFAULT_RESONANCE_TOLERANCE`.
pub fn is_body_resonant(body_id: i32, tolerance: f64) -> Result<bool> {
    let elements = minor_body_elements(body_id).ok_or(Error::UnknownBody { body_id })?;
    Ok(detect_mean_motion_resonance(&elements, tolerance).is_resonant)
}

/// Resonance information for a curated body using the default tolerance.
///
/// Unknown IDs return `None`. A known non-resonant body still returns `Some`
/// with `is_resonant=false`; inspect the result rather than treating presence
/// as a resonance verdict. This is a semi-major-axis ratio test, not proof of
/// libration of a resonant argument.
pub fn get_resonance_info(body_id: i32) -> Option<ResonanceResult> {
    let elements = minor_body_elements(body_id)?;
    Some(detect_mean_motion_resonance(
        &elements,
        DEFAULT_RESONANCE_TOLERANCE,
    ))
}

// ---------------------------------------------------------------------------
// Kepler solvers
// ---------------------------------------------------------------------------

/// Solve the elliptic Kepler equation `M = E - e*sin(E)`.
///
/// `m_rad` and the returned eccentric anomaly `E` are radians; eccentricity
/// `e` is dimensionless with physical domain `0 <= e < 1`. `tol` is the
/// absolute successive-iterate tolerance in radians. The finite input anomaly
/// takes a truncated remainder by one turn, then only positive values above
/// π fold downward; negative values below −π are retained.
///
/// Uses an `E = M` starter for `e < 0.8`, a Markley (1995) rational starter
/// otherwise, and at most 30 Newton steps. The last iterate is returned if
/// the tolerance is not reached. No domain or convergence error is reported;
/// invalid inputs, singular starters and overflow can yield nonfinite values.
pub fn solve_kepler_equation_elliptic(m_rad: f64, e: f64, tol: f64) -> f64 {
    let mut m = m_rad % (2.0 * std::f64::consts::PI);
    if m > std::f64::consts::PI {
        m -= 2.0 * std::f64::consts::PI;
    }
    let mut anomaly = if e < 0.8 {
        m
    } else {
        let alpha = (1.0 - e) / (4.0 * e + 0.5);
        let beta = m / (2.0 * (4.0 * e + 0.5));
        let aux = (beta * beta + alpha * alpha * alpha).sqrt();
        let z = (beta + aux).abs().powf(1.0 / 3.0);
        let s = z - alpha / z;
        let w = s - 0.078 * s * s * s * s * s / (1.0 + e);
        let e0 = m + e * (3.0 * w - 4.0 * w * w * w) / (1.0 + e * (1.0 - w * w * 1.5));
        if m < 0.0 { -e0.abs() } else { e0.abs() }
    };
    for _ in 0..30 {
        let f = anomaly - e * anomaly.sin() - m;
        let fp = 1.0 - e * anomaly.cos();
        let next = anomaly - f / fp;
        if (next - anomaly).abs() < tol {
            return next;
        }
        anomaly = next;
    }
    anomaly
}

/// Solve the hyperbolic Kepler equation `M = e*sinh(H) - H`.
///
/// `m_rad` is the hyperbolic mean-anomaly parameter; returns hyperbolic eccentric
/// anomaly `H`, conventionally in radians. Eccentricity is dimensionless with
/// physical domain `e > 1`; `tol` is the absolute successive-iterate tolerance.
/// At most 50 Newton steps run, stopping also when `|fp| < 1e-15`. The current
/// iterate is returned on either exhaustion or the derivative-floor stop,
/// without a convergence verdict. Inputs are not validated and nonfinite
/// arithmetic is not converted to an error.
pub fn solve_kepler_equation_hyperbolic(m_rad: f64, e: f64, tol: f64) -> f64 {
    let mut h = if e < 1.6 {
        if m_rad.abs() < std::f64::consts::PI {
            m_rad
        } else {
            m_rad.signum() * std::f64::consts::PI
        }
    } else if m_rad.abs() > 1e-10 {
        (2.0 * m_rad.abs() / e + 1.8).ln().copysign(m_rad)
    } else {
        0.0
    };
    for _ in 0..50 {
        let f = e * h.sinh() - h - m_rad;
        let fp = e * h.cosh() - 1.0;
        if fp.abs() < 1e-15 {
            break;
        }
        let next = h - f / fp;
        if (next - h).abs() < tol {
            h = next;
            break;
        }
        h = next;
    }
    h
}

/// Solve Barker's equation for parabolic orbits: closed-form `s = 2*sinh(asinh(3M)/3)`,
/// true anomaly `nu = 2*atan(s)`.
pub fn solve_barker_equation(m_rad: f64) -> f64 {
    let w = 3.0 * m_rad;
    let s = 2.0 * (w.asinh() / 3.0).sinh();
    2.0 * s.atan()
}

/// Select the Kepler solver for a dimensionless eccentricity `e`.
///
/// Within `|e - 1| < 1e-10`, Barker's equation returns true anomaly directly
/// in radians; otherwise the elliptic arm returns eccentric anomaly `E` and
/// the hyperbolic arm returns `H`. These are different anomaly definitions,
/// not interchangeable output angles. `tol` controls successive-iterate
/// convergence on the elliptic/hyperbolic arms and is unused by Barker's
/// closed form. The selected solver's unchecked-input and last-iterate
/// behavior applies; this dispatcher adds no validation or error result.
pub fn solve_kepler_equation(m_rad: f64, e: f64, tol: f64) -> f64 {
    if (e - 1.0).abs() < 1e-10 {
        solve_barker_equation(m_rad)
    } else if e < 1.0 {
        solve_kepler_equation_elliptic(m_rad, e, tol)
    } else {
        solve_kepler_equation_hyperbolic(m_rad, e, tol)
    }
}

// ---------------------------------------------------------------------------
// Closed-form corrections: Neptune longitude, resonant argument,
// libration, Yarkovsky, planet-element evolution
// ---------------------------------------------------------------------------

/// Neptune's mean longitude at `jd_tt`, degrees.
pub fn calc_neptune_mean_longitude(jd_tt: f64) -> f64 {
    (NEPTUNE_MEAN_LONGITUDE_J2000 + NEPTUNE_N * (jd_tt - J2000_EPOCH)).rem_euclid(360.0)
}

/// Resonant argument `phi = 3*lambda - 2*lambda_neptune - varpi` for a
/// plutino, degrees. `omega_pert`
/// and `m_pert` are the perturbed argument of perihelion and mean anomaly.
/// The longitude coefficients sum to zero, so `phi` is origin-invariant;
/// `phi` librates about 180 degrees for plutinos.
pub fn calc_resonant_argument_plutino(
    elements: &OrbitalElements,
    jd_tt: f64,
    omega_pert: f64,
    m_pert: f64,
) -> f64 {
    let lambda_tno = (elements.node + omega_pert + m_pert).rem_euclid(360.0);
    let varpi_tno = (elements.node + omega_pert).rem_euclid(360.0);
    let lambda_neptune = calc_neptune_mean_longitude(jd_tt);
    (3.0 * lambda_tno - 2.0 * lambda_neptune - varpi_tno).rem_euclid(360.0)
}

/// Longitude correction from resonant libration, degrees: `(amplitude / q) * sin(2*pi*dt / P +
/// phase)`, `0.0` for bodies without libration parameters.
pub fn calc_libration_correction(body_id: i32, jd_tt: f64) -> f64 {
    let Some(params) = libration_parameters(body_id) else {
        return 0.0;
    };
    let dt = jd_tt - J2000_EPOCH;
    let omega_lib = 2.0 * std::f64::consts::PI / params.period;
    (params.amplitude / (params.resonance_q as f64)) * (omega_lib * dt + params.phase_j2000).sin()
}

/// Yarkovsky semi-major-axis drift:
/// `a + (da/dt per day) * dt_days` for the six measured NEAs, `a`
/// otherwise. `da/dt` converts from AU/My with `1e6 * 365.25`.
pub fn apply_yarkovsky_correction(
    elements: &OrbitalElements,
    dt_days: f64,
    body_id: Option<i32>,
) -> f64 {
    let Some(id) = body_id else {
        return elements.a;
    };
    let Some(da_dt) = yarkovsky_da_dt_au_per_my(id) else {
        return elements.a;
    };
    elements.a + (da_dt / (1e6 * 365.25)) * dt_days
}

/// Giant-planet elements evolved to `jd_tt`: J2000.0 values plus the Simon et al.
/// linear century rates. Rows are
/// `(a, mass_ratio, e, omega_deg, node_deg, i_deg, neptune_threshold)`.
/// Time-evolved planets improve the forced eccentricity/inclination
/// vectors for century-scale propagations.
pub fn planet_elements_at_time(jd_tt: f64) -> [(f64, f64, f64, f64, f64, f64, f64); 4] {
    let t = (jd_tt - 2451545.0) / 36525.0;
    [
        (
            JUPITER_A,
            MASS_RATIO_JUPITER,
            JUPITER_E + JUPITER_DE * t,
            JUPITER_OMEGA + JUPITER_DOMEGA * t,
            JUPITER_NODE + JUPITER_DNODE * t,
            JUPITER_I + JUPITER_DI * t,
            0.0,
        ),
        (
            SATURN_A,
            MASS_RATIO_SATURN,
            SATURN_E + SATURN_DE * t,
            SATURN_OMEGA + SATURN_DOMEGA * t,
            SATURN_NODE + SATURN_DNODE * t,
            SATURN_I + SATURN_DI * t,
            0.0,
        ),
        (
            URANUS_A,
            MASS_RATIO_URANUS,
            URANUS_E + URANUS_DE * t,
            URANUS_OMEGA + URANUS_DOMEGA * t,
            URANUS_NODE + URANUS_DNODE * t,
            URANUS_I + URANUS_DI * t,
            0.0,
        ),
        (
            NEPTUNE_A,
            MASS_RATIO_NEPTUNE,
            NEPTUNE_E + NEPTUNE_DE * t,
            NEPTUNE_OMEGA + NEPTUNE_DOMEGA * t,
            NEPTUNE_NODE + NEPTUNE_DNODE * t,
            NEPTUNE_I + NEPTUNE_DI * t,
            20.0,
        ),
    ]
}

// ---------------------------------------------------------------------------
// Laplace-Lagrange secular core
// ---------------------------------------------------------------------------

/// Laplace coefficient `b_s^(j)(alpha)` by trapezoidal quadrature
/// : 500 steps above `alpha = 0.7`,
/// 200 below; `0.0` outside `(0,1)`; trapezoid endpoints weighted 1.0,
/// interior 2.0, prefactor `2/pi`. Method: Murray & Dermott §6.4;
/// Brouwer & Clemence Ch. XI.
pub fn laplace_coefficient(alpha: f64, s: f64, j: i32) -> f64 {
    if alpha >= 1.0 || alpha <= 0.0 {
        return 0.0;
    }
    let n_steps = if alpha > 0.7 { 500 } else { 200 };
    let dpsi = std::f64::consts::PI / (n_steps as f64);
    let mut result = 0.0;
    for k in 0..=n_steps {
        let psi = (k as f64) * dpsi;
        let denom = 1.0 - 2.0 * alpha * psi.cos() + alpha * alpha;
        if denom > 1e-10 {
            let integrand = ((j as f64) * psi).cos() / denom.powf(s);
            let weight = if k == 0 || k == n_steps { 1.0 } else { 2.0 };
            result += weight * integrand * dpsi / 2.0;
        }
    }
    2.0 * result / std::f64::consts::PI
}

/// Fallback secular rates `(d_omega, d_node, d_n)` in degrees/day: reduced Laplace-Lagrange
/// Jupiter/Saturn/Uranus/Neptune arms with the model interior/exterior
/// branch split (exact equality with a planet axis contributes nothing),
/// the project `d_n` apsidal-eccentricity coupling (`e < 0.99` guard) and
/// the bounded eight-item near-resonance amplification heuristic. See the
/// module header for the out-of-contract `e >= 1` division behavior.
/// Method: Murray & Dermott ch. 6-8; Brouwer & Clemence ch. XI, XVI.
/// Fallback only; SPK or ASSIST/REBOUND propagation takes scientific
/// precedence.
pub fn secular_perturbation_rates(elements: &OrbitalElements) -> (f64, f64, f64) {
    let a = elements.a;
    let e = elements.e;
    let i_rad = elements.i.to_radians();
    let n = elements.n;
    let n_rad = n.to_radians();
    let mut d_omega = 0.0;
    let mut d_node = 0.0;

    // Per-planet helper: (axis, mass_ratio) with the model branch rules.
    let arms = [
        (JUPITER_A, MASS_RATIO_JUPITER, true),
        (SATURN_A, MASS_RATIO_SATURN, true),
        (URANUS_A, MASS_RATIO_URANUS, true),
        (NEPTUNE_A, MASS_RATIO_NEPTUNE, false),
    ];
    for (axis, mu, always) in arms {
        // Neptune only for a > 20 AU (documented scope/performance choice).
        if !always && a <= 20.0 {
            continue;
        }
        // Interior arm (strict): unguarded 1/(1-e^2) like the model.
        if a < axis {
            let alpha = a / axis;
            let b32_1 = laplace_coefficient(alpha, 1.5, 1);
            let factor = n_rad * mu * alpha * b32_1;
            d_omega += (factor * (1.0 / 4.0) * (2.0 + 1.5 * e * e) / (1.0 - e * e)).to_degrees();
            d_node += (-factor * (1.0 / 4.0) * i_rad.cos() / (1.0 - e * e).sqrt()).to_degrees();
        } else if a > axis {
            // Exterior arm: max(0.01, 1-e^2) guards like the model.
            let alpha = axis / a;
            let b32_1 = laplace_coefficient(alpha, 1.5, 1);
            let factor = n_rad * mu * alpha * b32_1;
            d_omega +=
                (factor * (1.0 / 4.0) * (2.0 + 1.5 * e * e) / (1.0 - e * e).max(0.01)).to_degrees();
            d_node +=
                (-factor * (1.0 / 4.0) * i_rad.cos() / (1.0 - e * e).max(0.01).sqrt()).to_degrees();
        }
    }

    // Project mean-motion drift approximation (not a published equation):
    // d_n = degrees[-(3/2) * radians(d_omega) * e^2 / (1 - e^2)].
    let d_omega_rad = d_omega.to_radians();
    let d_n = if n_rad.abs() > 1e-20 && e < 0.99 {
        let e2 = e * e;
        (-1.5 * d_omega_rad * e2 / (1.0 - e2)).to_degrees()
    } else {
        0.0
    };

    // Project near-resonance validity heuristic: eight (n_planet, p, q,
    // max_amp) checks, dimensionless distance delta = |p*n - q*n_p| / n,
    // window 0.001 < delta < 0.05, scale min(1/(20*delta), max_amp).
    let resonance_checks = [
        (JUPITER_N, 3, 1, 3.0),
        (JUPITER_N, 5, 2, 2.5),
        (JUPITER_N, 7, 3, 2.0),
        (JUPITER_N, 2, 1, 3.0),
        (JUPITER_N, 3, 2, 3.0),
        (NEPTUNE_N, 2, 3, 3.0),
        (NEPTUNE_N, 1, 2, 2.5),
        (NEPTUNE_N, 3, 5, 2.0),
    ];
    let mut amplification: f64 = 1.0;
    for (n_planet, p_res, q_res, max_amp) in resonance_checks {
        if n.abs() > 1e-20 {
            let delta = ((p_res as f64) * n - (q_res as f64) * n_planet).abs() / n;
            if delta > 0.001 && delta < 0.05 {
                let amp = (1.0 / (delta * 20.0)).min(max_amp);
                amplification = amplification.max(amp);
            }
        }
    }
    if amplification > 1.0 {
        d_omega *= amplification;
        d_node *= amplification;
    }
    (d_omega, d_node, d_n)
}

/// Forced eccentricity/inclination vectors and proper secular frequencies.
///
/// Laplace-Lagrange theory splits a test particle's eccentricity vector
/// `(h,k) = (e sin(varpi), e cos(varpi))` and inclination vector
/// `(p,q) = (sin(i/2) sin(Omega), sin(i/2) cos(Omega))` into forced planetary
/// responses and free components rotating at frequencies `g` and `-s`.
/// Here `e` is eccentricity, `i` inclination, `Omega` ascending-node longitude,
/// and `varpi = omega + Omega` longitude of perihelion. Input element angles
/// are degrees, mean motion `n` is degrees/day, semi-major axis `a` is AU,
/// and `jd_tt` is a TT Julian Day. Planet elements evolve from J2000 using
/// the Simon et al. (1994) linear rates at that epoch.
///
/// For each planet `j`, let `mu_j` be its Sun-relative mass ratio,
/// `alpha = min(a/a_j, a_j/a)`, `D = max(0.01, 1-e^2)`, and let
/// `b_m = b_(3/2)^(m)(alpha)` be a Laplace coefficient. With `n` in radians/day,
/// accumulate `g_j = n mu_j alpha b_1 (2 + 1.5e^2)/(4D)` and
/// `s_j = n mu_j alpha b_1 cos(i)/(4 sqrt(D))`. The eccentricity coupling
/// `nu_j = n mu_j alpha b_2/4` multiplies the planet's `(h_j,k_j)` vector;
/// the inclination coupling `sigma_j = n mu_j alpha b_1/4` multiplies
/// `(p_j,q_j)`. Summing these numerators and dividing by `g = sum(g_j)` or
/// `s = sum(s_j)` gives the forced vectors. A denominator of magnitude at
/// most `1e-20` produces a zero forced vector rather than a singular ratio.
///
/// Planets excluded by the element table's interior threshold are skipped,
/// as are co-orbital separations inside one Hill radius
/// `a_j (mu_j/3)^(1/3)` and invalid `alpha` outside `(0,1)`. This secular
/// approximation does not resolve resonances or close encounters. For
/// `e > 0.3`, the eccentricity response is multiplied by
/// `1 + e^2/2 + 3e^4/8`; for `i > 20 degrees`, the inclination response
/// is multiplied by `1 + sin(i/2)^2`.
///
/// Returns `(g, h_forced, k_forced, s, p_forced, q_forced)`; frequencies are
/// radians/day and vector components are dimensionless. Inputs are assumed
/// to be finite physical orbital elements; this function does not validate
/// them or certify accuracy outside the secular approximation.
/// References: Murray & Dermott, *Solar System Dynamics*, sections 7.3–7.5,
/// equations 7.10 and 7.25–7.26; Brouwer & Clemence, *Methods of Celestial
/// Mechanics*, chapter XVI; Simon et al. (1994), A&A 282, 663.
pub fn forced_elements(elements: &OrbitalElements, jd_tt: f64) -> (f64, f64, f64, f64, f64, f64) {
    let a = elements.a;
    let e = elements.e;
    let i_rad = elements.i.to_radians();
    let n_rad = elements.n.to_radians();
    let e2 = e * e;
    let denom_e = (1.0 - e2).max(0.01);

    let mut g = 0.0;
    let mut s = 0.0;
    let mut h_num = 0.0;
    let mut k_num = 0.0;
    let mut p_num = 0.0;
    let mut q_num = 0.0;

    for (a_p, mu_p, e_p, omega_p_deg, node_p_deg, i_p_deg, a_thresh) in
        planet_elements_at_time(jd_tt)
    {
        if a_thresh > 0.0 && a < a_thresh {
            continue;
        }
        let hill_radius = a_p * (mu_p / 3.0).powf(1.0 / 3.0);
        if (a - a_p).abs() < hill_radius {
            continue;
        }
        let alpha = if a < a_p { a / a_p } else { a_p / a };
        if alpha >= 1.0 || alpha <= 0.0 {
            continue;
        }
        let b32_1 = laplace_coefficient(alpha, 1.5, 1);
        let b32_2 = laplace_coefficient(alpha, 1.5, 2);
        let base = n_rad * mu_p * alpha;

        g += base * b32_1 * (1.0 / 4.0) * (2.0 + 1.5 * e2) / denom_e;
        s += base * b32_1 * (1.0 / 4.0) * i_rad.cos() / denom_e.sqrt();
        let nu_j = base * b32_2 * (1.0 / 4.0);
        let sigma_j = base * b32_1 * (1.0 / 4.0);

        let varpi_p = (omega_p_deg + node_p_deg).to_radians();
        let i_p = i_p_deg.to_radians();
        let node_p = node_p_deg.to_radians();
        // Planet eccentricity/inclination vectors first, then scaled by the
        // coupling (`nu_j * h_j`, `sigma_j * p_j`): the model forms `h_j`,
        // `k_j`, `p_j`, `q_j` before the coupling multiply, and the fused
        // triple product rounds differently (1 ulp per planet, amplified by
        // the 1-second nod_aps stencil on near-circular orbits).
        let h_j = e_p * varpi_p.sin();
        let k_j = e_p * varpi_p.cos();
        let sin_half_ip = (i_p / 2.0).sin();
        let p_j = sin_half_ip * node_p.sin();
        let q_j = sin_half_ip * node_p.cos();
        h_num += nu_j * h_j;
        k_num += nu_j * k_j;
        p_num += sigma_j * p_j;
        q_num += sigma_j * q_j;
    }

    let (h_forced0, k_forced0) = if g.abs() > 1e-20 {
        (h_num / g, k_num / g)
    } else {
        (0.0, 0.0)
    };
    let (p_forced0, q_forced0) = if s.abs() > 1e-20 {
        (p_num / s, q_num / s)
    } else {
        (0.0, 0.0)
    };

    let (h_forced, k_forced) = if e > 0.3 {
        let corr = 1.0 + 0.5 * e2 + 0.375 * e2 * e2;
        (h_forced0 * corr, k_forced0 * corr)
    } else {
        (h_forced0, k_forced0)
    };
    let (p_forced, q_forced) = if elements.i > 20.0 {
        let sh = (i_rad / 2.0).sin();
        let corr = 1.0 + sh * sh;
        (p_forced0 * corr, q_forced0 * corr)
    } else {
        (p_forced0, q_forced0)
    };
    (g, h_forced, k_forced, s, p_forced, q_forced)
}

/// Perturbed elements at `jd_tt`:
/// `(omega_pert, node_pert, m_pert, n_pert, e_pert, i_pert)` with angles in
/// degrees. Short propagations (`|dt| < 1 d`) or `include_perturbations =
/// false` return the unperturbed propagation; `e_pert` clamps to
/// `[0, 0.999]` with no positive floor. First-order Laplace-Lagrange may
/// mislead for Neptune-resonant bodies such as plutinos and twotinos;
/// check `detect_mean_motion_resonance` / `is_body_resonant`.
pub fn apply_secular_perturbations(
    elements: &OrbitalElements,
    jd_tt: f64,
    include_perturbations: bool,
) -> (f64, f64, f64, f64, f64, f64) {
    let dt = jd_tt - elements.epoch;
    if !include_perturbations || dt.abs() < 1.0 {
        let m = (elements.m0 + elements.n * dt).rem_euclid(360.0);
        return (
            elements.omega,
            elements.node,
            m,
            elements.n,
            elements.e,
            elements.i,
        );
    }
    let (d_omega, d_node, d_n) = secular_perturbation_rates(elements);
    let omega_pert = (elements.omega + d_omega * dt).rem_euclid(360.0);
    let node_pert = (elements.node + d_node * dt).rem_euclid(360.0);
    let n_pert = elements.n + d_n;
    let m_pert = (elements.m0 + n_pert * dt).rem_euclid(360.0);

    let (g, h_forced, k_forced, s_freq, p_forced, q_forced) = forced_elements(elements, jd_tt);

    let varpi0 = (elements.omega + elements.node).to_radians();
    let h0 = elements.e * varpi0.sin();
    let k0 = elements.e * varpi0.cos();
    let dh = h0 - h_forced;
    let dk = k0 - k_forced;
    let e_free = (dh * dh + dk * dk).sqrt();
    let beta = if e_free > 1e-15 { dh.atan2(dk) } else { 0.0 };
    let h_t = h_forced + e_free * (g * dt + beta).sin();
    let k_t = k_forced + e_free * (g * dt + beta).cos();
    let e_pert = (h_t * h_t + k_t * k_t).sqrt().clamp(0.0, 0.999);

    let i0 = elements.i.to_radians();
    let node0 = elements.node.to_radians();
    let p0 = (i0 / 2.0).sin() * node0.sin();
    let q0 = (i0 / 2.0).sin() * node0.cos();
    let dp = p0 - p_forced;
    let dq = q0 - q_forced;
    let i_free = (dp * dp + dq * dq).sqrt();
    let gamma = if i_free > 1e-15 { dp.atan2(dq) } else { 0.0 };
    let p_t = p_forced + i_free * (-s_freq * dt + gamma).sin();
    let q_t = q_forced + i_free * (-s_freq * dt + gamma).cos();
    let sin_half_i_t = (p_t * p_t + q_t * q_t).sqrt().clamp(0.0, 1.0);
    let i_pert = (2.0 * sin_half_i_t.asin()).to_degrees();
    (omega_pert, node_pert, m_pert, n_pert, e_pert, i_pert)
}

// ---------------------------------------------------------------------------
// Short-period correction
// ---------------------------------------------------------------------------

/// First-order short-period Jupiter/Saturn correction, AU: `(dx, dy,
/// dz)` in the ecliptic J2000.0 frame. Angles are the perturbed
/// `(omega, node, i)` in radians. The dead `math.cos(...)` statement in
/// the model is omitted (no observable effect); the guarded
/// second-order term matches the shipped association. The dominant terms
/// oscillate at the synodic period. See also Murray & Dermott §6.9.
pub fn short_period_correction(
    elements: &OrbitalElements,
    jd_tt: f64,
    omega_rad: f64,
    node_rad: f64,
    i_rad: f64,
) -> (f64, f64, f64) {
    let a = elements.a;
    let e = elements.e;
    let n = elements.n;
    let dt = jd_tt - elements.epoch;
    if a <= 0.0 || e >= 1.0 {
        return (0.0, 0.0, 0.0);
    }
    let m_body = ((elements.m0 + n * dt).rem_euclid(360.0)).to_radians();
    let lambda_body = node_rad + omega_rad + m_body;

    // (axis, mass ratio, mean motion deg/day, mean longitude at J2000 rad).
    let sp_planets = [
        (
            JUPITER_A,
            MASS_RATIO_JUPITER,
            JUPITER_N,
            (JUPITER_NODE + JUPITER_OMEGA + 18.818).to_radians(),
        ),
        (
            SATURN_A,
            MASS_RATIO_SATURN,
            SATURN_N,
            (SATURN_NODE + SATURN_OMEGA + 320.346).to_radians(),
        ),
    ];
    let mut delta_lambda = 0.0;
    for (a_p, mu_p, n_p_deg, lambda_p_j2000) in sp_planets {
        let alpha = if a < a_p { a / a_p } else { a_p / a };
        if alpha >= 1.0 || alpha <= 0.01 {
            continue;
        }
        let n_p_rad = n_p_deg.to_radians();
        let lambda_p = lambda_p_j2000 + n_p_rad * (jd_tt - 2451545.0);
        let n_rad = n.to_radians();
        let dn = n_rad - n_p_rad;
        if dn.abs() < 1e-8 {
            continue;
        }
        let b12_1 = laplace_coefficient(alpha, 0.5, 1);
        let amplitude = 2.0 * mu_p * alpha * b12_1;
        delta_lambda += -amplitude * (lambda_body - lambda_p).sin() / (dn / n_rad);
        if (2.0 * n_rad - n_p_rad).abs() > 1e-8 {
            delta_lambda += -amplitude * e * (2.0 * m_body - lambda_p).sin()
                / ((2.0 * n_rad - n_p_rad) / n_rad)
                * 0.5;
        }
    }
    let r_approx = a * (1.0 - e * e);
    let (sin_lambda, cos_lambda) = lambda_body.sin_cos();
    let dx = -r_approx * sin_lambda * delta_lambda;
    let dy = r_approx * cos_lambda * delta_lambda * i_rad.cos();
    let dz = r_approx * cos_lambda * delta_lambda * i_rad.sin();
    (dx, dy, dz)
}

// ---------------------------------------------------------------------------
// Position core
// ---------------------------------------------------------------------------

/// Heliocentric ecliptic-J2000 position, AU: secular-perturbed elements for bound
/// orbits (with the differential plutino libration when `body_id` names a
/// modeled plutino), Kepler solve per orbit class (parabolic rows carry
/// `q` in `a`), Yarkovsky drift for elliptic orbits, perifocal rotation
/// with the perturbed angles, short-period add-on for perturbed bound
/// orbits. Expected agreement with the model is belt ~10-30 arcsec,
/// TNO ~1-3 arcmin, and modeled plutinos within 2 deg/century; a
/// non-finite `jd_tt` is [`Error::NonFiniteArgument`]. Method:
/// Murray & Dermott ch. 7-8; Malhotra (1995); Curtis §3; Vallado §2.3.
pub fn calc_minor_body_position(
    elements: &OrbitalElements,
    jd_tt: f64,
    include_perturbations: bool,
    body_id: Option<i32>,
) -> Result<[f64; 3]> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let e = elements.e;
    let dt = jd_tt - elements.epoch;

    let (omega_pert, node_pert, m_deg, e_pert, i_pert) = if include_perturbations && e < 1.0 {
        let (om, no, m, _, ee, ii) = apply_secular_perturbations(elements, jd_tt, true);
        (om, no, m, ee, ii)
    } else {
        (
            elements.omega,
            elements.node,
            (elements.m0 + elements.n * dt).rem_euclid(360.0),
            e,
            elements.i,
        )
    };
    // For hyperbolic rows the perturbed eccentricity is unused.
    let e_use = if e < 1.0 { e_pert } else { e };

    let m_rad: f64 = if (e - 1.0).abs() < 1e-10 {
        (elements.m0 + elements.n * dt).to_radians()
    } else if e < 1.0 {
        let mut m_corr = m_deg;
        if let Some(id) = body_id {
            if libration_parameters(id).is_some() {
                m_corr += calc_libration_correction(id, jd_tt)
                    - calc_libration_correction(id, elements.epoch);
            }
        }
        m_corr.rem_euclid(360.0).to_radians()
    } else {
        (elements.m0 + elements.n * dt).to_radians()
    };
    let anomaly = solve_kepler_equation(m_rad, e_use, 1e-8);

    let (nu, r) = if (e - 1.0).abs() < 1e-10 {
        let p = 2.0 * elements.a;
        (anomaly, p / (1.0 + anomaly.cos()))
    } else if e < 1.0 {
        let nu = 2.0
            * ((1.0 + e_use).sqrt() * (anomaly / 2.0).sin())
                .atan2((1.0 - e_use).sqrt() * (anomaly / 2.0).cos());
        let a_use = apply_yarkovsky_correction(elements, dt, body_id);
        (nu, a_use * (1.0 - e_use * anomaly.cos()))
    } else {
        let nu = 2.0
            * ((e + 1.0).sqrt() * (anomaly / 2.0).sinh())
                .atan2((e - 1.0).sqrt() * (anomaly / 2.0).cosh());
        (nu, elements.a.abs() * (e * anomaly.cosh() - 1.0))
    };

    let (x_orb, y_orb) = (r * nu.cos(), r * nu.sin());
    let (omega_rad, node_rad) = if e < 1.0 {
        (omega_pert.to_radians(), node_pert.to_radians())
    } else {
        (elements.omega.to_radians(), elements.node.to_radians())
    };
    let i_rad = i_pert.to_radians();
    let (sin_om, cos_om) = omega_rad.sin_cos();
    let (sin_no, cos_no) = node_rad.sin_cos();
    let (sin_i, cos_i) = i_rad.sin_cos();
    let p11 = cos_om * cos_no - sin_om * sin_no * cos_i;
    let p12 = -sin_om * cos_no - cos_om * sin_no * cos_i;
    let p21 = cos_om * sin_no + sin_om * cos_no * cos_i;
    let p22 = -sin_om * sin_no + cos_om * cos_no * cos_i;
    let p31 = sin_om * sin_i;
    let p32 = cos_om * sin_i;
    let mut x = p11 * x_orb + p12 * y_orb;
    let mut y = p21 * x_orb + p22 * y_orb;
    let mut z = p31 * x_orb + p32 * y_orb;

    if include_perturbations && e < 1.0 {
        let (dx, dy, dz) = short_period_correction(elements, jd_tt, omega_rad, node_rad, i_rad);
        x += dx;
        y += dy;
        z += dz;
    }
    Ok([x, y, z])
}

// ---------------------------------------------------------------------------
// Table lookup: base + multi-epoch nodes with midpoint cross-fade
// ---------------------------------------------------------------------------

use crate::constants::{AST_OFFSET, PHOLUS};
use crate::minor_body_tables::{MINOR_BODY_MULTI_TABLE, MINOR_BODY_TABLE, MinorBodyRow};
use std::collections::HashMap;

/// Half-width of the cross-fade window around epoch midpoints, days
pub const EPOCH_BLEND_HALF_DAYS: f64 = 182.625;
/// Epochs closer than this are the same node when weaving the curated
/// base entry into the multi-epoch table.
pub const EPOCH_MERGE_TOL_DAYS: f64 = 1.0;

fn elements_from_row(row: &MinorBodyRow) -> OrbitalElements {
    OrbitalElements {
        name: row.name.to_string(),
        epoch: row.epoch,
        a: row.a,
        e: row.e,
        i: row.i_deg,
        omega: row.omega_deg,
        node: row.node_deg,
        m0: row.m0_deg,
        n: row.n_deg_day,
        physical_spice_id: crate::constants::get_spk_body_info_from_map(row.body_id)
            .map(|(_, target)| target),
    }
}

/// Curated single-epoch elements for a body id: `None` for unknown ids.
pub fn minor_body_elements(body_id: i32) -> Option<OrbitalElements> {
    MINOR_BODY_TABLE
        .iter()
        .find(|row| row.body_id == body_id)
        .map(elements_from_row)
}

/// Curated element rows in their defining table insertion order.
/// Each result owns its scalar fields and retains its row's private physical
/// identity. No identifier is inferred from names. Units and epoch basis are the same as
/// [`minor_body_elements`]. This performs no propagation or acquisition.
pub fn ordered_minor_body_elements() -> Vec<(i32, OrbitalElements)> {
    MINOR_BODY_TABLE
        .iter()
        .map(|row| (row.body_id, elements_from_row(row)))
        .collect()
}

/// Multi-epoch element rows in their original defining table order.
/// The base table remains distinct: this does not merge, deduplicate or sort
/// epochs. Each result owns its elements with the same units and epoch basis
/// as [`minor_body_elements`]. No model evaluation or source acquisition runs.
pub fn ordered_minor_body_multi_elements() -> Vec<(i32, OrbitalElements)> {
    MINOR_BODY_MULTI_TABLE
        .iter()
        .map(|row| (row.body_id, elements_from_row(row)))
        .collect()
}

/// Base + multi-epoch nodes for a body in ascending epoch order: the curated base entry is woven in unless a
/// table node already sits at its epoch. `None` without a multi-epoch
/// table. Recomputed per call (tens of rows); the model caches the
/// immutable merge, which is observationally identical. The curated
/// epoch is the closest present-decade fit: dropping it loses ~300
/// arcsec near the present, and cross-fading avoids steps at node
/// switches.
pub fn merged_epoch_nodes(body_id: i32) -> Option<Vec<OrbitalElements>> {
    let base = minor_body_elements(body_id)?;
    let mut nodes: Vec<OrbitalElements> = MINOR_BODY_MULTI_TABLE
        .iter()
        .filter(|row| row.body_id == body_id)
        .map(elements_from_row)
        .collect();
    if nodes.is_empty() {
        return None;
    }
    if nodes
        .iter()
        .all(|e| (base.epoch - e.epoch).abs() > EPOCH_MERGE_TOL_DAYS)
    {
        nodes.push(base);
    }
    nodes.sort_by(|a, b| {
        a.epoch
            .partial_cmp(&b.epoch)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    Some(nodes)
}

/// Nearest element epoch plus a cross-fade partner near midpoints
/// : `(primary, secondary, weight)`
/// with `weight` in `(0, 0.5]`, else `(primary, None, 0.0)`. Unknown ids
/// are [`Error::UnknownBody`]. First-in-list wins exact midpoint ties,
/// clipping to the physical upper endpoint.
pub fn epoch_elements_blend(
    body_id: i32,
    jd_tt: f64,
) -> Result<(OrbitalElements, Option<OrbitalElements>, f64)> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let base = minor_body_elements(body_id).ok_or(Error::UnknownBody { body_id })?;
    let Some(entries) = merged_epoch_nodes(body_id) else {
        return Ok((base, None, 0.0));
    };
    let mut nearest_idx = 0;
    let mut best = f64::INFINITY;
    for (idx, e) in entries.iter().enumerate() {
        let dist = (jd_tt - e.epoch).abs();
        if dist < best {
            best = dist;
            nearest_idx = idx;
        }
    }
    let nearest = &entries[nearest_idx];
    let other = if jd_tt >= nearest.epoch && nearest_idx + 1 < entries.len() {
        Some(&entries[nearest_idx + 1])
    } else if jd_tt < nearest.epoch && nearest_idx > 0 {
        Some(&entries[nearest_idx - 1])
    } else {
        None
    };
    let Some(other) = other else {
        return Ok((entries[nearest_idx].clone(), None, 0.0));
    };
    let midpoint = (nearest.epoch + other.epoch) / 2.0;
    let dist = (jd_tt - midpoint).abs();
    if dist >= EPOCH_BLEND_HALF_DAYS {
        return Ok((entries[nearest_idx].clone(), None, 0.0));
    }
    let weight = 0.5 * (1.0 - dist / EPOCH_BLEND_HALF_DAYS);
    Ok((entries[nearest_idx].clone(), Some(other.clone()), weight))
}

/// Best elements for a body at the target time: the blend primary. The
/// nearest-node choice caps two-body propagation at ~5 years between
/// 10-year nodes.
pub fn closest_epoch_elements(body_id: i32, jd_tt: f64) -> Result<OrbitalElements> {
    Ok(epoch_elements_blend(body_id, jd_tt)?.0)
}

// ---------------------------------------------------------------------------
// Heliocentric entries
// ---------------------------------------------------------------------------

/// Cartesian to spherical `(lon_deg, lat_deg, dist)`:
/// unclamped `asin(z / r)`; a non-finite quotient is
/// [`Error::NonFiniteArgument`].
pub fn cartesian_to_spherical(x: f64, y: f64, z: f64) -> Result<(f64, f64, f64)> {
    let r = (x * x + y * y + z * z).sqrt();
    let q = z / r;
    if !q.is_finite() {
        return Err(Error::NonFiniteArgument { value: q });
    }
    Ok((
        y.atan2(x).to_degrees().rem_euclid(360.0),
        q.asin().to_degrees(),
        r,
    ))
}

/// Range-type provider misses that fall through to the Keplerian path
fn is_fallthrough_miss(err: &Error) -> bool {
    matches!(
        err,
        Error::EphemerisRange { .. }
            | Error::SpkOutOfRange { .. }
            | Error::LebOutOfRange { .. }
            | Error::UnknownBody { .. }
            | Error::InvalidBody { .. }
    )
}

fn spk_or_none(
    spk: Option<SpkPositionProvider<'_>>,
    body_id: i32,
    jd_tt: f64,
) -> Result<Option<(f64, f64, f64)>> {
    let Some(lookup) = spk else {
        return Ok(None);
    };
    match lookup(body_id, jd_tt) {
        Ok(hit) => Ok(hit),
        Err(err) if is_fallthrough_miss(&err) => Ok(None),
        Err(err) => Err(err),
    }
}

/// Heliocentric ecliptic coordinates, degrees/AU: SPK position when the
/// provider serves the body, else the epoch-blended Keplerian propagation
/// (cross-faded in Cartesian space) shaped to spherical, precessed from
/// J2000 to date, plus `dpsi_rad` for the true-of-date longitude. Unknown
/// ids are [`Error::UnknownBody`]. Served SPK
/// bodies agree at the sub-arcsecond level; the Keplerian fallback is
/// ~10-30 arcsec in the belt and ~1-3 arcmin for TNOs, growing
/// off-epoch. Heliocentric only: geocentric callers subtract Earth's
/// heliocentric position.
pub fn calc_minor_body_heliocentric(
    body_id: i32,
    jd_tt: f64,
    use_spk: bool,
    spk_position: Option<SpkPositionProvider<'_>>,
    dpsi_rad: f64,
) -> Result<(f64, f64, f64)> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    if !dpsi_rad.is_finite() {
        return Err(Error::NonFiniteArgument { value: dpsi_rad });
    }
    calc_minor_body_heliocentric_with_nutation(body_id, jd_tt, use_spk, spk_position, &|| {
        Ok(dpsi_rad)
    })
}

/// Heliocentric minor-body place with nutation acquired at its final consumed
/// stage. Coordinates and blending follow [`calc_minor_body_heliocentric`].
/// Unknown bodies refuse before SPK access; a served SPK consumes no nutation.
/// The fallback shapes and precesses the Cartesian blend before invoking the
/// callback for radians of nutation longitude. Existing propagation, frame and
/// provider errors remain terminal before that callback. TT JD is checked by
/// the caller/propagator; no eager nutation admission is added here.
pub fn calc_minor_body_heliocentric_with_nutation(
    body_id: i32,
    jd_tt: f64,
    use_spk: bool,
    spk_position: Option<SpkPositionProvider<'_>>,
    nutation: &dyn Fn() -> Result<f64>,
) -> Result<(f64, f64, f64)> {
    if minor_body_elements(body_id).is_none() {
        return Err(Error::UnknownBody { body_id });
    }
    if use_spk {
        if let Some(hit) = spk_or_none(spk_position, body_id, jd_tt)? {
            return Ok(hit);
        }
    }
    let (primary, secondary, weight) = epoch_elements_blend(body_id, jd_tt)?;
    let [x, y, z] = calc_minor_body_position(&primary, jd_tt, true, Some(body_id))?;
    let (x, y, z) = match secondary {
        Some(other) if weight > 0.0 => {
            let [xb, yb, zb] = calc_minor_body_position(&other, jd_tt, true, Some(body_id))?;
            (
                x * (1.0 - weight) + xb * weight,
                y * (1.0 - weight) + yb * weight,
                z * (1.0 - weight) + zb * weight,
            )
        }
        _ => (x, y, z),
    };
    shape_heliocentric_with_nutation(x, y, z, jd_tt, nutation)
}

/// Spherical shaping shared by the heliocentric entries: Cartesian to
/// spherical, J2000-to-date precession, nutation-in-longitude add-on.
fn shape_heliocentric(
    x: f64,
    y: f64,
    z: f64,
    jd_tt: f64,
    dpsi_rad: f64,
) -> Result<(f64, f64, f64)> {
    shape_heliocentric_with_nutation(x, y, z, jd_tt, &|| Ok(dpsi_rad))
}

/// Shape an acquired J2000 ecliptic Cartesian vector, in AU, to an of-date
/// longitude/latitude/radius in degrees/degrees/AU. Nutation radians are
/// consumed only after the existing spherical and precession operations.
/// Their checked domains/errors remain intact; no model or provider is read.
pub fn shape_heliocentric_with_nutation(
    x: f64,
    y: f64,
    z: f64,
    jd_tt: f64,
    nutation: &dyn Fn() -> Result<f64>,
) -> Result<(f64, f64, f64)> {
    let (lon, lat, r) = cartesian_to_spherical(x, y, z)?;
    let (lon_p, lat_p) = crate::hypothetical::precess_ecliptic(lon, lat, J2000_EPOCH, jd_tt)?;
    Ok((
        (lon_p + nutation()?.to_degrees()).rem_euclid(360.0),
        lat_p,
        r,
    ))
}

/// Bodies with dedicated ids shadowing their catalog numbers.
pub const SPECIAL_ASTEROID_IDS: [(i64, i32); 6] = [
    (2060, CHIRON),
    (5145, PHOLUS),
    (1, CERES),
    (2, PALLAS),
    (3, JUNO),
    (4, VESTA),
];

/// Resolution of an asteroid catalog number without network: `Ok(Some(id))`
/// for predefined-table and special-mapped numbers, `Ok(None)` when only
/// the SBDB leg could serve it. Non-positive numbers are
/// [`Error::InvalidBody`].
pub fn resolve_asteroid_number(asteroid_number: i64) -> Result<Option<i32>> {
    if asteroid_number <= 0 {
        return Err(Error::InvalidBody {
            body_id: i32::try_from(asteroid_number).unwrap_or(i32::MIN),
            operation: "calc_asteroid_by_number",
        });
    }
    if let Ok(as_i32) = i32::try_from(asteroid_number.saturating_add(AST_OFFSET as i64)) {
        if minor_body_elements(as_i32).is_some() {
            return Ok(Some(as_i32));
        }
    }
    for (number, mapped) in SPECIAL_ASTEROID_IDS {
        if number == asteroid_number {
            return Ok(Some(mapped));
        }
    }
    Ok(None)
}

/// Position for any numbered asteroid:
/// predefined table, then the six special ids, then the SPK provider when
/// `use_spk`, then the SBDB provider. An unserved number is
/// [`Error::UnknownBody`]; sealed refusals propagate while other provider
/// failures fall back to the unknown-body outcome exactly like the model
/// (whose fetch swallows transport errors into `None`). The SBDB leg
/// propagates with `body_id = None`, so no libration correction applies.
/// SPK-served numbers agree at the sub-arcsecond level; the Keplerian
/// path is ~10-30 arcsec in the belt and ~1-3 arcmin for TNOs. SBDB
/// elements arrive through the injected provider; element caching lives
/// in the SBDB layer.
#[allow(clippy::too_many_arguments)]
pub fn calc_asteroid_by_number(
    asteroid_number: i64,
    jd_tt: f64,
    use_spk: bool,
    spk_position: Option<SpkPositionProvider<'_>>,
    sbdb_elements: Option<SbdbElementsProvider<'_>>,
    dpsi_rad: f64,
) -> Result<(f64, f64, f64)> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    if !dpsi_rad.is_finite() {
        return Err(Error::NonFiniteArgument { value: dpsi_rad });
    }
    if let Some(body_id) = resolve_asteroid_number(asteroid_number)? {
        return calc_minor_body_heliocentric(body_id, jd_tt, use_spk, spk_position, dpsi_rad);
    }
    let unknown = Error::UnknownBody {
        body_id: i32::try_from(asteroid_number.saturating_add(AST_OFFSET as i64))
            .unwrap_or(i32::MAX),
    };
    if use_spk {
        if let Ok(as_i32) = i32::try_from(asteroid_number.saturating_add(AST_OFFSET as i64)) {
            if let Some(hit) = spk_or_none(spk_position, as_i32, jd_tt)? {
                return Ok(hit);
            }
        }
    }
    let Some(fetch) = sbdb_elements else {
        return Err(unknown);
    };
    match fetch(asteroid_number) {
        Ok(Some(elements)) => {
            let [x, y, z] = calc_minor_body_position(&elements, jd_tt, true, None)?;
            shape_heliocentric(x, y, z, jd_tt, dpsi_rad)
        }
        Ok(None) => Err(unknown),
        Err(err) if matches!(err, Error::NetworkSealed { .. }) => Err(err),
        Err(_) => Err(unknown),
    }
}

// ---------------------------------------------------------------------------
// JPL SBDB lookup layer (request shape + narrow response parsers)
// ---------------------------------------------------------------------------

/// JPL Small-Body Database API endpoint.
pub const SBDB_API_URL: &str = "https://ssd-api.jpl.nasa.gov/sbdb.api";
/// the model `User-Agent` for SBDB requests.
pub const SBDB_USER_AGENT: &str = "libephemeris/1.0 (Python)";

/// Elements-request URL:
/// `sstr=<number>&phys-par=false&full-prec=true`.
pub fn sbdb_elements_url(asteroid_number: i64) -> String {
    format!("{SBDB_API_URL}?sstr={asteroid_number}&phys-par=false&full-prec=true")
}

/// Name-lookup URL:
/// `sstr=<quoted name>&phys-par=false`, `urllib.quote` default (`safe='/'`).
pub fn sbdb_name_url(name: &str) -> String {
    format!(
        "{SBDB_API_URL}?sstr={}&phys-par=false",
        percent_encode_query(name)
    )
}

/// Percent-encode a query value like `urllib.parse.quote` with the default
/// `safe='/'`: unreserved `A-Z a-z 0-9 _ . - ~` plus `/` pass through,
/// every other UTF-8 byte becomes `%XX` (uppercase hex).
pub fn percent_encode_query(value: &str) -> String {
    let mut out = String::new();
    for byte in value.bytes() {
        match byte {
            b'A'..=b'Z' | b'a'..=b'z' | b'0'..=b'9' | b'_' | b'.' | b'-' | b'~' | b'/' => {
                out.push(byte as char)
            }
            _ => out.push_str(&format!("%{byte:02X}")),
        }
    }
    out
}

/// Minimal JSON value for the narrow SBDB shapes (objects, arrays,
/// strings with escapes, numbers kept raw, literals). Depth-guarded.
#[derive(Clone, Debug, PartialEq)]
enum JsonValue {
    Null,
    Bool(bool),
    Num(String),
    Str(String),
    Arr(Vec<JsonValue>),
    Obj(Vec<(String, JsonValue)>),
}

impl JsonValue {
    fn get(&self, key: &str) -> Option<&JsonValue> {
        if let JsonValue::Obj(entries) = self {
            entries.iter().find(|(k, _)| k == key).map(|(_, v)| v)
        } else {
            None
        }
    }
    /// A number, numeric string or Boolean SBDB field.
    fn as_f64(&self) -> Option<f64> {
        match self {
            JsonValue::Num(raw) => raw.parse::<f64>().ok(),
            JsonValue::Str(text) => text.parse::<f64>().ok(),
            JsonValue::Bool(true) => Some(1.0),
            JsonValue::Bool(false) => Some(0.0),
            _ => None,
        }
    }
    fn as_str(&self) -> Option<&str> {
        if let JsonValue::Str(text) = self {
            Some(text)
        } else {
            None
        }
    }
    /// Raw token for SPK-ID decoding (number or string).
    fn as_raw(&self) -> Option<&str> {
        match self {
            JsonValue::Num(raw) | JsonValue::Str(raw) => Some(raw),
            _ => None,
        }
    }
}

struct JsonParser<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> JsonParser<'a> {
    fn new(text: &'a str) -> Self {
        Self {
            bytes: text.as_bytes(),
            pos: 0,
        }
    }
    fn ws(&mut self) {
        while self.pos < self.bytes.len() && self.bytes[self.pos].is_ascii_whitespace() {
            self.pos += 1;
        }
    }
    fn peek(&self) -> Option<u8> {
        self.bytes.get(self.pos).copied()
    }
    fn eat(&mut self, byte: u8) -> bool {
        if self.peek() == Some(byte) {
            self.pos += 1;
            true
        } else {
            false
        }
    }
    fn parse_value(&mut self, depth: usize) -> Option<JsonValue> {
        if depth > 64 {
            return None;
        }
        self.ws();
        match self.peek()? {
            b'{' => self.parse_object(depth),
            b'[' => self.parse_array(depth),
            b'"' => Some(JsonValue::Str(self.parse_string()?)),
            b't' => self.parse_literal("true", JsonValue::Bool(true)),
            b'f' => self.parse_literal("false", JsonValue::Bool(false)),
            b'n' => self.parse_literal("null", JsonValue::Null),
            _ => self.parse_number(),
        }
    }
    fn parse_literal(&mut self, word: &str, value: JsonValue) -> Option<JsonValue> {
        if self.bytes[self.pos..].starts_with(word.as_bytes()) {
            self.pos += word.len();
            Some(value)
        } else {
            None
        }
    }
    fn parse_object(&mut self, depth: usize) -> Option<JsonValue> {
        self.eat(b'{');
        let mut entries = Vec::new();
        self.ws();
        if self.eat(b'}') {
            return Some(JsonValue::Obj(entries));
        }
        loop {
            self.ws();
            if self.peek() != Some(b'"') {
                return None;
            }
            let key = self.parse_string()?;
            self.ws();
            if !self.eat(b':') {
                return None;
            }
            let value = self.parse_value(depth + 1)?;
            entries.push((key, value));
            self.ws();
            if self.eat(b',') {
                continue;
            }
            if self.eat(b'}') {
                return Some(JsonValue::Obj(entries));
            }
            return None;
        }
    }
    fn parse_array(&mut self, depth: usize) -> Option<JsonValue> {
        self.eat(b'[');
        let mut items = Vec::new();
        self.ws();
        if self.eat(b']') {
            return Some(JsonValue::Arr(items));
        }
        loop {
            items.push(self.parse_value(depth + 1)?);
            self.ws();
            if self.eat(b',') {
                continue;
            }
            if self.eat(b']') {
                return Some(JsonValue::Arr(items));
            }
            return None;
        }
    }
    fn hex4(&mut self) -> Option<u32> {
        if self.pos + 4 > self.bytes.len() {
            return None;
        }
        let mut code: u32 = 0;
        for i in 0..4 {
            code = code * 16 + (self.bytes[self.pos + i] as char).to_digit(16)?;
        }
        self.pos += 4;
        Some(code)
    }
    fn parse_string(&mut self) -> Option<String> {
        if !self.eat(b'"') {
            return None;
        }
        let mut out = String::new();
        loop {
            let byte = *self.bytes.get(self.pos)?;
            self.pos += 1;
            match byte {
                b'"' => return Some(out),
                b'\\' => {
                    let esc = *self.bytes.get(self.pos)?;
                    self.pos += 1;
                    match esc {
                        b'"' => out.push('"'),
                        b'\\' => out.push('\\'),
                        b'/' => out.push('/'),
                        b'b' => out.push('\u{0008}'),
                        b'f' => out.push('\u{000C}'),
                        b'n' => out.push('\n'),
                        b'r' => out.push('\r'),
                        b't' => out.push('\t'),
                        b'u' => {
                            let mut code = self.hex4()?;
                            if (0xD800..0xDC00).contains(&code) {
                                if self.bytes.get(self.pos) == Some(&b'\\')
                                    && self.bytes.get(self.pos + 1) == Some(&b'u')
                                {
                                    self.pos += 2;
                                    let low = self.hex4()?;
                                    if (0xDC00..0xE000).contains(&low) {
                                        code = 0x10000 + ((code - 0xD800) << 10) + (low - 0xDC00);
                                    } else {
                                        return None;
                                    }
                                } else {
                                    return None;
                                }
                            }
                            out.push(char::from_u32(code)?);
                        }
                        _ => return None,
                    }
                }
                _ => {
                    // Raw UTF-8 passthrough (JSON source is UTF-8).
                    let rest = std::str::from_utf8(&self.bytes[self.pos - 1..]).ok()?;
                    let ch = rest.chars().next()?;
                    if (ch as u32) < 0x20 {
                        return None;
                    }
                    out.push(ch);
                    self.pos += ch.len_utf8() - 1;
                }
            }
        }
    }
    fn parse_number(&mut self) -> Option<JsonValue> {
        let start = self.pos;
        if self.peek() == Some(b'-') {
            self.pos += 1;
        }
        match self.peek()? {
            b'0' => self.pos += 1,
            b'1'..=b'9' => {
                while matches!(self.peek(), Some(b'0'..=b'9')) {
                    self.pos += 1;
                }
            }
            _ => return None,
        }
        if self.peek() == Some(b'.') {
            self.pos += 1;
            if !matches!(self.peek(), Some(b'0'..=b'9')) {
                return None;
            }
            while matches!(self.peek(), Some(b'0'..=b'9')) {
                self.pos += 1;
            }
        }
        if matches!(self.peek(), Some(b'e' | b'E')) {
            self.pos += 1;
            if matches!(self.peek(), Some(b'+' | b'-')) {
                self.pos += 1;
            }
            if !matches!(self.peek(), Some(b'0'..=b'9')) {
                return None;
            }
            while matches!(self.peek(), Some(b'0'..=b'9')) {
                self.pos += 1;
            }
        }
        let raw = std::str::from_utf8(&self.bytes[start..self.pos]).ok()?;
        Some(JsonValue::Num(raw.to_string()))
    }
}

fn parse_sbdb_document(text: &str) -> Option<JsonValue> {
    let mut parser = JsonParser::new(text);
    let value = parser.parse_value(0)?;
    parser.ws();
    if parser.pos != parser.bytes.len() {
        return None;
    }
    Some(value)
}

/// Decode an SBDB `spkid` token to an asteroid catalog number: modern `20000000+N` first, then legacy
/// `2000000+N`, else the bare id; non-positive results are `None`.
/// Decimal-string ids are additionally accepted by truncation within an
/// absolute bound; the model integer conversion returns `None` there.
pub fn decode_sbdb_spkid(raw: &str) -> Option<i64> {
    let id = raw.parse::<i64>().ok().or_else(|| {
        raw.parse::<f64>().ok().and_then(|v| {
            if v.is_finite() && v.abs() < 1e15 {
                Some(v.trunc() as i64)
            } else {
                None
            }
        })
    })?;
    let number = if id >= 20_000_000 {
        id - 20_000_000
    } else if id >= 2_000_000 {
        id - 2_000_000
    } else {
        id
    };
    if number <= 0 { None } else { Some(number) }
}

/// Retain explicit SBDB identity without changing numerical-response admission.
/// A response identifier must be an exact integer and agree with a positive
/// request selector. Missing response identity can retain that selector;
/// contradictory or malformed identity stays anonymous, never name-inferred.
fn sbdb_physical_spice_id(document: &JsonValue, asteroid_number: i64) -> Option<i32> {
    let requested = (asteroid_number > 0).then_some(asteroid_number);
    let number = match document
        .get("object")
        .and_then(|object| object.get("spkid"))
    {
        Some(value) => {
            let raw = value.as_raw()?;
            // The general lookup decoder permits fractional-string truncation.
            // Such a token cannot establish a physical self-force identity.
            raw.parse::<i64>().ok()?;
            let response = decode_sbdb_spkid(raw)?;
            if requested.is_some_and(|request| request != response) {
                return None;
            }
            response
        }
        None => requested?,
    };
    // Selected force files use the physical legacy NAIF asteroid identity,
    // independently of modern SBDB's twenty-million-offset representation.
    let target = number.checked_add(i64::from(crate::constants::NAIF_ASTEROID_OFFSET))?;
    i32::try_from(target).ok()
}

/// Parse an SBDB elements response into [`OrbitalElements`]: top-level
/// `"error"` or an empty element list is `None`; the epoch must be
/// nonzero; `a`/`e`/`i`/`n` are required while `w`/`om`/`ma` default to
/// `0.0`; the name falls back to `"Asteroid <number>"`. Non-finite
/// values are treated as missing. An unambiguous explicit catalogue identity
/// is retained privately for N-body self-force selection; malformed or
/// conflicting identity metadata does not change numerical parsing.
pub fn parse_sbdb_elements_response(text: &str, asteroid_number: i64) -> Option<OrbitalElements> {
    let doc = parse_sbdb_document(text)?;
    let top = if let JsonValue::Obj(_) = doc {
        &doc
    } else {
        return None;
    };
    if top.get("error").is_some() {
        return None;
    }
    let orbit = top.get("orbit")?;
    let elements = if let JsonValue::Arr(items) = orbit.get("elements")? {
        items
    } else {
        return None;
    };
    if elements.is_empty() {
        return None;
    }
    let mut map: HashMap<&str, f64> = HashMap::new();
    for entry in elements {
        let name = entry.get("name")?.as_str()?;
        let value = entry.get("value")?.as_f64()?;
        map.insert(name, value);
    }
    let epoch = orbit.get("epoch")?.as_f64()?;
    if epoch == 0.0 {
        return None;
    }
    let (a, e, i, n) = (
        *map.get("a")?,
        *map.get("e")?,
        *map.get("i")?,
        *map.get("n")?,
    );
    let name = top
        .get("object")
        .and_then(|obj| obj.get("fullname"))
        .and_then(JsonValue::as_str)
        .map(str::to_string)
        .unwrap_or_else(|| format!("Asteroid {asteroid_number}"));
    let mut elements = OrbitalElements::new(
        name,
        epoch,
        a,
        e,
        i,
        *map.get("w").unwrap_or(&0.0),
        *map.get("om").unwrap_or(&0.0),
        *map.get("ma").unwrap_or(&0.0),
        n,
    )
    .ok()?;
    elements.physical_spice_id = sbdb_physical_spice_id(top, asteroid_number);
    Some(elements)
}

/// Parse an SBDB name-lookup response into a catalog number: `spkid`, then all-digit `des`,
/// then the leading token of `fullname`; non-positive results are `None`.
pub fn parse_sbdb_name_response(text: &str) -> Option<i64> {
    let doc = parse_sbdb_document(text)?;
    let top = if let JsonValue::Obj(_) = doc {
        &doc
    } else {
        return None;
    };
    if top.get("error").is_some() {
        return None;
    }
    let obj = top.get("object")?;
    if let Some(raw) = obj.get("spkid").and_then(JsonValue::as_raw) {
        return decode_sbdb_spkid(raw);
    }
    if let Some(des) = obj.get("des").and_then(JsonValue::as_str) {
        if !des.is_empty() && des.bytes().all(|b| b.is_ascii_digit()) {
            let number: i64 = des.parse().ok()?;
            return if number > 0 { Some(number) } else { None };
        }
    }
    let fullname = obj
        .get("fullname")
        .and_then(JsonValue::as_str)
        .unwrap_or("");
    let first = fullname.split_whitespace().next().unwrap_or("");
    if !first.is_empty() && first.bytes().all(|b| b.is_ascii_digit()) {
        let number: i64 = first.parse().ok()?;
        return if number > 0 { Some(number) } else { None };
    }
    None
}

/// Caller-owned map of asteroid catalogue numbers to retained SBDB elements.
///
/// This cache performs no fetching, expiry or size-based eviction. Successful
/// fetches can populate it; callers choose its lifetime and when to clear it.
/// `get` borrows a stored row, `insert` takes ownership and replaces any prior
/// row, and cloning the cache makes an independent copy rather than a shared
/// mutable resource. Stored epochs and units follow [`OrbitalElements`].
#[derive(Clone, Debug, Default)]
pub struct SbdbElementsCache {
    entries: HashMap<i64, OrbitalElements>,
}

impl SbdbElementsCache {
    /// Empty caller-owned cache.
    pub fn new() -> Self {
        Self::default()
    }
    /// Cached elements for `asteroid_number`, if present.
    pub fn get(&self, asteroid_number: i64) -> Option<&OrbitalElements> {
        self.entries.get(&asteroid_number)
    }
    /// Store `elements` under `asteroid_number`, replacing any prior row.
    pub fn insert(&mut self, asteroid_number: i64, elements: OrbitalElements) {
        self.entries.insert(asteroid_number, elements);
    }
    /// Number of cached rows.
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    /// Whether the cache holds no rows.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
    /// Clear the cache, returning the cleared entry count.
    pub fn clear(&mut self) -> usize {
        let count = self.entries.len();
        self.entries.clear();
        count
    }
}

/// Well-known asteroid names and numbers: lowercase name to catalog number. "Europa" and
/// "Pandora" name the main-belt asteroids, not the moons; "lilith" names
/// main-belt (1181) Lilith, not the lunar apogee point.
pub const LOCAL_ASTEROID_NAMES: &[(&str, i64)] = &[
    ("chiron", 2060),
    ("pholus", 5145),
    ("ceres", 1),
    ("pallas", 2),
    ("juno", 3),
    ("vesta", 4),
    ("eris", 136199),
    ("sedna", 90377),
    ("haumea", 136108),
    ("makemake", 136472),
    ("ixion", 28978),
    ("orcus", 90482),
    ("quaoar", 50000),
    ("nessus", 7066),
    ("asbolus", 8405),
    ("chariklo", 10199),
    ("gonggong", 225088),
    ("varuna", 20000),
    ("apophis", 99942),
    ("flora", 8),
    ("eunomia", 15),
    ("hygiea", 10),
    ("interamnia", 704),
    ("davida", 511),
    ("europa", 52),
    ("sylvia", 87),
    ("psyche", 16),
    ("eros", 433),
    ("amor", 1221),
    ("icarus", 1566),
    ("toro", 1685),
    ("sappho", 80),
    ("pandora", 55),
    ("lilith", 1181),
    ("hidalgo", 944),
    ("toutatis", 4179),
    ("itokawa", 25143),
    ("bennu", 101955),
    ("ryugu", 162173),
];

/// Caller-owned map of normalized asteroid names to catalogue numbers.
///
/// Keys use Unicode trimming and lowercase conversion. Local defaults are
/// seeded only while the map is empty; inserting a custom entry first prevents
/// that automatic seeding. `get_local` can therefore mutate an empty cache.
/// Remote successes added by [`get_asteroid_number`] have no expiry or automatic
/// eviction. Clones own independent maps; clearing permits the next local
/// lookup to seed defaults again. No lookup here performs network acquisition.
#[derive(Clone, Debug, Default)]
pub struct AsteroidNameCache {
    entries: HashMap<String, i64>,
}

impl AsteroidNameCache {
    /// Empty caller-owned cache (seed with [`AsteroidNameCache::seed_local`]).
    pub fn new() -> Self {
        Self::default()
    }
    /// Seed the well-known local mappings; a no-op once non-empty
    pub fn seed_local(&mut self) {
        if !self.entries.is_empty() {
            return;
        }
        for (name, number) in LOCAL_ASTEROID_NAMES {
            self.entries.insert((*name).to_string(), *number);
        }
    }
    /// Look up a trimmed, lowercased name, lazily seeding an empty map.
    ///
    /// Empty normalized names return `None` without seeding. A nonempty map
    /// is never supplemented with missing default rows by this operation.
    pub fn get_local(&mut self, name: &str) -> Option<i64> {
        let key = name.trim().to_lowercase();
        if key.is_empty() {
            return None;
        }
        self.seed_local();
        self.entries.get(&key).copied()
    }
    /// Store `asteroid_number` under the model-normalized name form
    /// (`trim.to_lowercase`), replacing any prior row.
    pub fn insert(&mut self, name: &str, asteroid_number: i64) {
        // the model caches under the normalized form (`name_lower`).
        self.entries
            .insert(name.trim().to_lowercase(), asteroid_number);
    }
    /// Number of cached name mappings.
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    /// Whether the cache holds no mappings.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
    /// Clear the cache, returning the cleared entry count.
    pub fn clear(&mut self) -> usize {
        let count = self.entries.len();
        self.entries.clear();
        count
    }
}

/// Name-to-number lookup: local mappings
/// first, then the remote provider. Remote successes are cached under the
/// normalized name; anything else is `None`. Sealed refusals propagate.
/// Timeout, TLS and request-purpose handling belong to the injected remote
/// provider, which receives the raw (unnormalized) name; the SBDB name-URL
/// helper is available for building such a request.
pub fn get_asteroid_number(
    name: &str,
    cache: &mut AsteroidNameCache,
    remote: Option<AsteroidNameProvider<'_>>,
) -> Result<Option<i64>> {
    let key = name.trim().to_lowercase();
    if key.is_empty() {
        return Ok(None);
    }
    cache.seed_local();
    if let Some(number) = cache.entries.get(&key).copied() {
        return Ok(Some(number));
    }
    let Some(fetch) = remote else {
        return Ok(None);
    };
    match fetch(name) {
        Ok(Some(number)) => {
            cache.entries.insert(key, number);
            Ok(Some(number))
        }
        Ok(None) => Ok(None),
        Err(err) if matches!(err, Error::NetworkSealed { .. }) => Err(err),
        Err(_) => Ok(None),
    }
}

// ---------------------------------------------------------------------------
// Horizons-SPK provisioning policy (pure half)
// ---------------------------------------------------------------------------

/// Earliest Horizons SPK-generation date.
pub const HORIZONS_SPK_DATE_MIN: &str = "1600-01-01";
/// Latest Horizons SPK-generation date.
pub const HORIZONS_SPK_DATE_MAX: &str = "2500-01-01";
/// `HORIZONS_SPK_DATE_MIN` as a Julian Day (UTC).
pub const HORIZONS_SPK_JD_MIN: f64 = 2305448.5;
/// `HORIZONS_SPK_DATE_MAX` as a Julian Day (UTC).
pub const HORIZONS_SPK_JD_MAX: f64 = 2634167.5;

/// Clamp a requested SPK window into the Horizons generation range
/// : `None` when the
/// request lies entirely outside the padded range; non-finite or
/// negative padding is [`Error::InvalidPadding`]. Non-finite window
/// endpoints are [`Error::NonFiniteArgument`]; partial overlaps clamp to
/// the padded range.
pub fn clamp_spk_window(
    jd_start: f64,
    jd_end: f64,
    boundary_padding_days: f64,
) -> Result<Option<(f64, f64)>> {
    if !jd_start.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_start });
    }
    if !jd_end.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_end });
    }
    if !boundary_padding_days.is_finite() || boundary_padding_days < 0.0 {
        return Err(Error::InvalidPadding {
            entry: "auto_download_asteroid_spk",
            value: boundary_padding_days,
        });
    }
    let clamp_start = HORIZONS_SPK_JD_MIN - boundary_padding_days;
    let clamp_end = HORIZONS_SPK_JD_MAX + boundary_padding_days;
    if jd_start < clamp_start || jd_end > clamp_end {
        let start = jd_start.max(clamp_start);
        let end = jd_end.min(clamp_end);
        if start >= end {
            return Ok(None);
        }
        return Ok(Some((start, end)));
    }
    Ok(Some((jd_start, jd_end)))
}

/// SPK download target for a body: the major-asteroid row first, then
/// `SPK_BODY_NAME_MAP` with the `spk_body_name` display name (or the bare
/// id); `None` for unsupported bodies. The live download/registration
/// stays with the SPK/NET families. Major-row Horizons ids use the `Name;`
/// form because bare numerics collide with planet-barycenter ids; the map
/// leg carries zero in the catalog-number slot.
pub fn spk_download_target(body_id: i32) -> Option<(i32, &'static str, i32, String)> {
    if let Some((number, horizons_id, naif_id, name)) = crate::spk::get_major_asteroid_info(body_id)
    {
        return Some((number, horizons_id, naif_id, name.to_string()));
    }
    if let Some((horizons_id, naif_id)) = crate::constants::get_spk_body_info_from_map(body_id) {
        let name = crate::spk::spk_body_name(body_id)
            .map(str::to_string)
            .unwrap_or_else(|| format!("{body_id}"));
        return Some((0, horizons_id, naif_id, name));
    }
    None
}

/// Major asteroids with SPK support as `(body_id, display name)`.
pub fn list_major_asteroids() -> Vec<(i32, &'static str)> {
    crate::spk::MAJOR_ASTEROID_SPK_INFO
        .iter()
        .map(|(id, _, _, _, name)| (*id, *name))
        .collect()
}

/// Whether an SPK kernel is registered for a body over a caller registry predicate.
pub fn is_spk_available_for_body(is_registered: &dyn Fn(i32) -> bool, body_id: i32) -> bool {
    is_registered(body_id)
}

// ---------------------------------------------------------------------------
// SBDB fetch + Horizons auto-download compositions (API-07 U 16 g; the model
// `minor_bodies.fetch_orbital_elements_from_sbdb`,
// `auto_download_asteroid_spk` and `ensure_major_asteroid_spk` over
// explicit handles, no process-global state).
// ---------------------------------------------------------------------------

/// Fetch SBDB elements with the caller-owned cache: a cache hit answers without
/// socket work; a non-positive number is `None` without socket work
/// ; the sealed policy is
/// [`Error::NetworkSealed`] before any socket work; every other
/// transport/parse failure is `None`.
/// Successes populate `cache` (failures never do). `timeout_secs`
/// feeds the single `timeout` into both transport budgets
/// (invalid values are [`Error::InvalidTimeout`]).
pub fn fetch_sbdb_elements(
    asteroid_number: i64,
    timeout_secs: f64,
    policy: crate::context::EffectiveNetworkPolicy,
    cache: &mut SbdbElementsCache,
) -> Result<Option<OrbitalElements>> {
    if let Some(hit) = cache.get(asteroid_number) {
        return Ok(Some(hit.clone()));
    }
    if asteroid_number <= 0 {
        return Ok(None);
    }
    let config = crate::net::TransportConfig::with_timeout_secs(timeout_secs)?;
    crate::context::require_network(
        policy,
        &format!("SBDB orbital elements for asteroid {asteroid_number}"),
    )?;
    let url = sbdb_elements_url(asteroid_number);
    let request = crate::net::HttpRequest::get(&url)?.header("User-Agent", SBDB_USER_AGENT);
    let purpose = format!("SBDB orbital elements for asteroid {asteroid_number}");
    let response = match crate::net::fetch(&request, &purpose, policy, &config) {
        Ok(response) => response,
        Err(Error::NetworkSealed { .. }) => {
            return Err(Error::NetworkSealed { purpose });
        }
        Err(_) => return Ok(None),
    };
    let elements = response
        .body_text()
        .ok()
        .and_then(|text| parse_sbdb_elements_response(text, asteroid_number));
    match elements {
        Some(elements) => {
            cache.insert(asteroid_number, elements.clone());
            Ok(Some(elements))
        }
        None => Ok(None),
    }
}

/// Legacy wall-clock Julian-day estimate; despite its name, this is not TT.
///
/// Converts elapsed Unix seconds with a fixed 86400-second day and the Unix
/// epoch offset. No leap-second table, TAI offset or TT conversion is applied,
/// so the result is a Unix/UTC-like calendar estimate, not a precision time
/// scale. If the system clock precedes the Unix epoch, the duration failure
/// substitutes zero elapsed seconds and returns the Unix epoch Julian day.
/// Callers needing TT must use an explicit time-scale conversion instead.
pub fn current_jd_tt() -> f64 {
    const UNIX_TO_JD: f64 = 2_440_587.5;
    UNIX_TO_JD
        + std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_secs_f64() / 86_400.0)
            .unwrap_or(0.0)
}

/// TT Julian Day as a `YYYY-MM-DD` UTC calendar date: TT→TAI by the exact offset, TAI→UTC
/// over the packaged leap table, then the proleptic-Gregorian
/// reconstruction. `None` on any conversion failure (the caller then
/// reports the unserved window, never an invented date).
pub fn tt_jd_to_utc_date(jd_tt: f64) -> Option<String> {
    let tai = crate::time::tt_to_tai_jd(jd_tt).ok()?;
    let table = crate::time::hardcoded_leap_seconds();
    let utc = crate::time::tai_to_utc_jd(tai, &table).ok()?;
    let date = crate::time::revjul(utc, crate::time::Calendar::Gregorian).ok()?;
    if !(1..=12).contains(&date.month) || !(1..=31).contains(&date.day) {
        return None;
    }
    Some(format!(
        "{:04}-{:02}-{:02}",
        date.year, date.month, date.day
    ))
}

/// Horizons SPK NAIF id for a downloaded kernel: the catalog number matched against the kernel's own targets
/// first (exact 20000000+N vs 2000000+N convention), then the lone
/// small-body target, then the lone non-center target, then the
/// legacy name deduction. `None` means the kernel cannot serve `ipl`
pub fn deduce_downloaded_naif_id(
    targets: &[i32],
    body: &str,
    asteroid_number: Option<i64>,
) -> Option<i32> {
    if let Some(number) = asteroid_number {
        if let Some(hit) = crate::spk::find_naif_id_for_asteroid(targets, number) {
            return Some(hit);
        }
    }
    if let Some(hit) = crate::spk::detect_naif_id_from_targets(targets, None) {
        return Some(hit);
    }
    let unique: Vec<i32> = {
        let mut seen = targets
            .iter()
            .copied()
            .filter(|t| *t != 0 && *t != 10)
            .collect::<Vec<_>>();
        seen.sort_unstable();
        seen.dedup();
        seen
    };
    if let [single] = unique.as_slice() {
        return Some(*single);
    }
    crate::spk::deduce_naif_id(body, None)
}

/// Automatically download an SPK kernel for a minor body:
/// unsupported bodies and entirely out-of-window requests are
/// `Ok(None)`; an already-registered body answers its stored path
/// (membership only, like the check — never a
/// coverage probe); otherwise the window (explicit, else now ± 10 y)
/// clamps via [`clamp_spk_window`], renders to UTC dates, downloads
/// through [`crate::spk::download_horizons_spk_to`] into `dir` under
/// the generated cache filename, and registers the kernel on the
/// file-deduced NAIF id. The sealed policy is [`Error::NetworkSealed`]
/// before any socket work; every other download/registration failure
/// is `Ok(None)`. `timeout_secs`
/// feeds both transport budgets ([`Error::InvalidTimeout`] when
/// invalid).
#[allow(clippy::too_many_arguments)]
pub fn auto_download_asteroid_spk_to(
    body_id: i32,
    jd_start: Option<f64>,
    jd_end: Option<f64>,
    force: bool,
    boundary_padding_days: f64,
    dir: &std::path::Path,
    registry: &mut crate::spk::SpkRegistry,
    policy: crate::context::EffectiveNetworkPolicy,
    timeout_secs: f64,
) -> Result<Option<std::path::PathBuf>> {
    let Some((_number, horizons_id, _naif_id, name)) = spk_download_target(body_id) else {
        return Ok(None);
    };
    if !force {
        if let Some(info) = registry.body_info(body_id) {
            return Ok(Some(std::path::PathBuf::from(info.kernel.clone())));
        }
    }
    let now = current_jd_tt();
    let start = jd_start.unwrap_or(now - 3652.5);
    let end = jd_end.unwrap_or(now + 3652.5);
    let Some((from, to)) = clamp_spk_window(start, end, boundary_padding_days)? else {
        return Ok(None);
    };
    let config = crate::net::TransportConfig::with_timeout_secs(timeout_secs)?;
    let (Some(start_date), Some(end_date)) = (tt_jd_to_utc_date(from), tt_jd_to_utc_date(to))
    else {
        return Ok(None);
    };
    let filename = crate::spk::generate_spk_cache_filename(horizons_id, from, to);
    let dest = dir.join(filename);
    let outcome = (|| {
        crate::spk::download_horizons_spk_to(
            horizons_id,
            &start_date,
            &end_date,
            "500@0",
            &dest,
            true,
            policy,
            &config,
        )?;
        let key = registry.register_file(&dest)?;
        let targets = crate::spk::spk_targets_of(registry, &key)?;
        let number = crate::spk::extract_asteroid_number(horizons_id);
        let naif = deduce_downloaded_naif_id(&targets, horizons_id, number)
            .ok_or(Error::SpkBodyNotFound { target: body_id })?;
        registry.register_body(body_id, &key, naif)?;
        Ok(dest.clone())
    })();
    match outcome {
        Ok(path) => Ok(Some(path)),
        Err(Error::NetworkSealed { .. }) => Err(Error::NetworkSealed {
            purpose: format!("Auto SPK download for {name}"),
        }),
        Err(_) => Ok(None),
    }
}

/// Ensure an SPK is available for a minor body:
/// `true` when the registry already serves the body, else one
/// [`auto_download_asteroid_spk_to`] attempt centered on `jd`
/// (`None` centers on now) with zero padding. Sealed refusal
/// propagates; download failures answer `false`.
pub fn ensure_major_asteroid_spk_to(
    body_id: i32,
    jd: Option<f64>,
    dir: &std::path::Path,
    registry: &mut crate::spk::SpkRegistry,
    policy: crate::context::EffectiveNetworkPolicy,
    timeout_secs: f64,
) -> Result<bool> {
    if registry.body_info(body_id).is_some() {
        return Ok(true);
    }
    let (start, end) = match jd {
        Some(center) => {
            if !center.is_finite() {
                return Err(Error::NonFiniteArgument { value: center });
            }
            (Some(center - 3652.5), Some(center + 3652.5))
        }
        None => (None, None),
    };
    Ok(auto_download_asteroid_spk_to(
        body_id,
        start,
        end,
        false,
        0.0,
        dir,
        registry,
        policy,
        timeout_secs,
    )?
    .is_some())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{
        APOPHIS as C_APOPHIS, BENNU as C_BENNU, CERES as C_CERES, CHIRON as C_CHIRON,
        ERIS as C_ERIS, IXION as C_IXION, JUNO as C_JUNO, PALLAS as C_PALLAS, VESTA as C_VESTA,
    };

    fn sample_elements() -> OrbitalElements {
        OrbitalElements::new(
            "Test".to_string(),
            2461000.5,
            2.77,
            0.08,
            10.5,
            73.0,
            80.0,
            231.0,
            0.214,
        )
        .expect("finite sample")
    }

    #[test]
    fn catalogue_identity_survives_cloning_but_scalar_records_never_infer_it() {
        let catalogue = minor_body_elements(C_CERES).expect("curated Ceres");
        assert_eq!(catalogue.physical_spice_id(), Some(2_000_001));
        assert_eq!(catalogue.clone().physical_spice_id(), Some(2_000_001));
        let copy = OrbitalElements::new(
            catalogue.name().into(),
            catalogue.epoch(),
            catalogue.a(),
            catalogue.e(),
            catalogue.i(),
            catalogue.omega(),
            catalogue.node(),
            catalogue.m0(),
            catalogue.n(),
        )
        .unwrap();
        assert_eq!(copy.physical_spice_id(), None);
        let source_record = OrbitalElements::from_nbody_source_record(
            catalogue.name().into(),
            catalogue.epoch(),
            catalogue.a(),
            catalogue.e(),
            catalogue.i(),
            catalogue.omega(),
            catalogue.node(),
            catalogue.m0(),
            catalogue.n(),
        );
        assert_eq!(source_record.physical_spice_id(), None);
        let ordered = ordered_minor_body_elements();
        let retained = ordered.iter().find(|(body, _)| *body == C_CERES).unwrap();
        assert_eq!(retained.1.physical_spice_id(), Some(2_000_001));
        assert!(
            merged_epoch_nodes(C_CERES)
                .unwrap()
                .iter()
                .all(|row| row.physical_spice_id() == Some(2_000_001))
        );
    }

    #[test]
    fn base_table_dimensions_and_spot_rows() {
        assert_eq!(MINOR_BODY_TABLE.len(), 37);
        let ceres = minor_body_elements(C_CERES).expect("ceres");
        assert_eq!(ceres.name(), "Ceres");
        assert_eq!(ceres.epoch(), 2461000.5);
        assert_eq!(ceres.a(), 2.765615651508659);
        assert_eq!(ceres.e(), 0.07957631994408416);
        assert_eq!(ceres.i(), 10.58788658206854);
        assert_eq!(ceres.omega(), 73.29975464616518);
        assert_eq!(ceres.node(), 80.24963090816965);
        assert_eq!(ceres.m0(), 231.5397330043706);
        assert_eq!(ceres.n(), 0.2142971214271186);
        let sedna = minor_body_elements(90377 + AST_OFFSET).expect("sedna");
        assert_eq!(sedna.name(), "Sedna");
        assert_eq!(sedna.a(), 549.5405279241455);
        assert_eq!(sedna.e(), 0.8612974378772624);
        assert_eq!(sedna.n(), 7.650758332945998e-05);
        assert!(minor_body_elements(999999).is_none());
    }

    #[test]
    fn multi_table_dimensions_and_order() {
        assert_eq!(MINOR_BODY_MULTI_TABLE.len(), 2916);
        // 36 bodies x 81 nodes; Bennu has no multi-epoch table.
        assert!(merged_epoch_nodes(C_BENNU).is_none());
        let nodes = merged_epoch_nodes(C_CERES).expect("ceres nodes");
        // Base woven in: 81 table nodes + 1 base node.
        assert_eq!(nodes.len(), 82);
        for pair in nodes.windows(2) {
            assert!(pair[0].epoch < pair[1].epoch);
        }
        let first = &MINOR_BODY_MULTI_TABLE[0];
        assert_eq!(first.name, "Amor");
        assert_eq!(first.epoch, 2323707.5);
        assert_eq!(first.a, 1.92613152020203);
        assert_eq!(first.e, 0.4357364028758751);
        assert_eq!(first.i_deg, 12.09373176563);
        assert_eq!(first.omega_deg, 20.45683431043);
        assert_eq!(first.node_deg, 174.491443655);
        assert_eq!(first.m0_deg, 120.9384564683);
        assert_eq!(first.n_deg_day, 0.36870172239378984);
    }

    #[test]
    fn blend_midpoint_and_far_field() {
        // Far from any midpoint: no partner.
        let (primary, secondary, weight) =
            epoch_elements_blend(C_CERES, 2461000.5).expect("blend at base epoch");
        assert_eq!(primary.epoch(), 2461000.5);
        assert!(secondary.is_none());
        assert_eq!(weight, 0.0);
        // Midpoint between the first two Ceres nodes gets weight 0.5.
        let nodes = merged_epoch_nodes(C_CERES).expect("nodes");
        let midpoint = (nodes[0].epoch + nodes[1].epoch) / 2.0;
        let (primary, secondary, weight) =
            epoch_elements_blend(C_CERES, midpoint).expect("midpoint blend");
        assert_eq!(primary.epoch(), nodes[0].epoch);
        assert_eq!(secondary.expect("partner").epoch(), nodes[1].epoch);
        assert_eq!(weight, 0.5);
        // Bennu: base only.
        let (primary, secondary, weight) = epoch_elements_blend(C_BENNU, 2461000.5).expect("bennu");
        assert_eq!(primary.name(), "Bennu");
        assert!(secondary.is_none());
        assert_eq!(weight, 0.0);
        // Unknown id refuses.
        assert!(matches!(
            epoch_elements_blend(999999, 2461000.5),
            Err(Error::UnknownBody { .. })
        ));
    }

    #[test]
    fn kepler_analytic_pins() {
        // M = 0 is a fixed point for every class.
        assert_eq!(solve_kepler_equation(0.0, 0.0, 1e-8), 0.0);
        assert_eq!(solve_kepler_equation(0.0, 0.9, 1e-8), 0.0);
        assert_eq!(solve_barker_equation(0.0), 0.0);
        assert_eq!(solve_kepler_equation_hyperbolic(0.0, 1.5, 1e-8), 0.0);
        // Circular: E = M.
        let e = solve_kepler_equation(1.234, 0.0, 1e-8);
        assert!((e - 1.234).abs() < 1e-12);
        // Elliptic residual at a hard case.
        let e_star = solve_kepler_equation(2.0, 0.95, 1e-8);
        assert!((e_star - 0.95 * e_star.sin() - 2.0).abs() < 1e-7);
        // Hyperbolic residual.
        let h = solve_kepler_equation(1.0, 2.0, 1e-8);
        assert!((2.0 * h.sinh() - h - 1.0).abs() < 1e-7);
        // Dispatcher routes the parabolic class to Barker.
        assert_eq!(
            solve_kepler_equation(0.7, 1.0, 1e-8),
            solve_barker_equation(0.7)
        );
    }

    #[test]
    fn resonance_detection_shapes() {
        let ixion = minor_body_elements(C_IXION).expect("ixion");
        let result = detect_mean_motion_resonance(&ixion, DEFAULT_RESONANCE_TOLERANCE);
        assert!(result.is_resonant);
        let info = result.resonance.expect("plutino info");
        assert_eq!((info.p, info.q), (2, 3));
        assert_eq!(info.name, "plutino");
        assert!(result.warning_message.expect("warning").contains("plutino"));
        let ceres = minor_body_elements(C_CERES).expect("ceres");
        let quiet = detect_mean_motion_resonance(&ceres, DEFAULT_RESONANCE_TOLERANCE);
        assert!(!quiet.is_resonant);
        assert!(quiet.resonance.is_none());
        assert!(quiet.warning_message.is_none());
        assert!(is_body_resonant(C_IXION, DEFAULT_RESONANCE_TOLERANCE).expect("ixion"));
        assert!(!is_body_resonant(C_CERES, DEFAULT_RESONANCE_TOLERANCE).expect("ceres"));
        assert!(matches!(
            is_body_resonant(999999, DEFAULT_RESONANCE_TOLERANCE),
            Err(Error::UnknownBody { .. })
        ));
        assert!(get_resonance_info(999999).is_none());
    }

    #[test]
    fn libration_and_yarkovsky_spots() {
        assert_eq!(calc_libration_correction(777777, 2461000.5), 0.0);
        assert!(!has_libration_model(777777));
        assert!(has_libration_model(C_IXION));
        let params = get_libration_parameters(C_IXION).expect("ixion params");
        assert_eq!((params.resonance_p, params.resonance_q), (2, 3));
        let elements = sample_elements();
        assert_eq!(
            apply_yarkovsky_correction(&elements, 100.0, None),
            elements.a
        );
        assert_eq!(
            apply_yarkovsky_correction(&elements, 100.0, Some(777777)),
            elements.a
        );
        // Apophis drifts inward: da/dt negative.
        let drifted = apply_yarkovsky_correction(&elements, 36525.0, Some(C_APOPHIS));
        assert!(drifted < elements.a);
    }

    #[test]
    fn neptune_longitude_j2000() {
        assert_eq!(
            calc_neptune_mean_longitude(J2000_EPOCH),
            NEPTUNE_MEAN_LONGITUDE_J2000
        );
    }

    #[test]
    fn validate_warning_shapes() {
        let bad_a = OrbitalElements::new(
            "Bad".to_string(),
            2461000.5,
            -1.0,
            0.1,
            10.0,
            0.0,
            0.0,
            0.0,
            0.1,
        )
        .expect("finite");
        let warnings = validate_elements(&bad_a);
        assert_eq!(warnings.len(), 1);
        assert!(warnings[0].contains("non-positive"));
        let bad_e = OrbitalElements::new(
            "Bad".to_string(),
            2461000.5,
            1.0,
            1.5,
            200.0,
            0.0,
            0.0,
            0.0,
            0.1,
        )
        .expect("finite");
        assert_eq!(validate_elements(&bad_e).len(), 2);
        assert!(validate_elements(&sample_elements()).is_empty());
        assert!(
            OrbitalElements::new(
                "Bad".to_string(),
                f64::NAN,
                1.0,
                0.1,
                10.0,
                0.0,
                0.0,
                0.0,
                0.1
            )
            .is_err()
        );
    }

    #[test]
    fn spherical_shaping_pins() {
        assert_eq!(
            cartesian_to_spherical(1.0, 0.0, 0.0).expect("x-axis"),
            (0.0, 0.0, 1.0)
        );
        assert_eq!(
            cartesian_to_spherical(0.0, 1.0, 0.0).expect("y-axis"),
            (90.0, 0.0, 1.0)
        );
        assert!(cartesian_to_spherical(0.0, 0.0, 0.0).is_err());
    }

    #[test]
    fn percent_encode_vectors() {
        assert_eq!(percent_encode_query("Ceres"), "Ceres");
        assert_eq!(percent_encode_query("2060"), "2060");
        assert_eq!(percent_encode_query("a b/c+d"), "a%20b/c%2Bd");
        assert_eq!(percent_encode_query("Éris"), "%C3%89ris");
        assert!(sbdb_elements_url(433).ends_with("?sstr=433&phys-par=false&full-prec=true"));
        assert!(sbdb_name_url("Apophis").ends_with("?sstr=Apophis&phys-par=false"));
    }

    #[test]
    fn sbdb_elements_parser_shapes() {
        let payload = r#"{"object": {"fullname": "433 Eros (1898 DQ)"}, "orbit": {"epoch": 2461000.5, "elements": [{"name": "e", "value": 0.2229}, {"name": "a", "value": "1.458"}, {"name": "q", "value": 1.133}, {"name": "i", "value": 10.83}, {"name": "om", "value": 304.3}, {"name": "w", "value": 178.8}, {"name": "ma", "value": 10.0}, {"name": "n", "value": 0.559}]}}"#;
        let elements = parse_sbdb_elements_response(payload, 433).expect("eros");
        assert_eq!(elements.name(), "433 Eros (1898 DQ)");
        assert_eq!(elements.a(), 1.458);
        assert_eq!(elements.node(), 304.3);
        // Error key refuses.
        assert!(parse_sbdb_elements_response(r#"{"error": "not found"}"#, 1).is_none());
        // Empty element list refuses.
        assert!(
            parse_sbdb_elements_response(r#"{"orbit": {"epoch": 2461000.5, "elements": []}}"#, 1)
                .is_none()
        );
        // Zero epoch refuses.
        assert!(
            parse_sbdb_elements_response(
                r#"{"orbit": {"epoch": 0, "elements": [{"name": "a", "value": 1.0}]}}"#,
                1
            )
            .is_none()
        );
        // Missing core element refuses.
        assert!(
            parse_sbdb_elements_response(
                r#"{"orbit": {"epoch": 2461000.5, "elements": [{"name": "a", "value": 1.0}]}}"#,
                1
            )
            .is_none()
        );
        // Defaults for w/om/ma + fallback name.
        let minimal = parse_sbdb_elements_response(
            r#"{"orbit": {"epoch": 2461000.5, "elements": [{"name": "a", "value": 2.0}, {"name": "e", "value": 0.1}, {"name": "i", "value": 5.0}, {"name": "n", "value": 0.3}]}}"#,
            12345,
        )
        .expect("minimal");
        assert_eq!(minimal.name(), "Asteroid 12345");
        assert_eq!(
            (minimal.omega(), minimal.node(), minimal.m0()),
            (0.0, 0.0, 0.0)
        );
        // Malformed JSON refuses.
        assert!(parse_sbdb_elements_response("{oops", 1).is_none());
        assert!(parse_sbdb_elements_response("[1, 2]", 1).is_none());
    }

    #[test]
    fn sbdb_elements_parser_retains_only_unambiguous_catalogue_identity() {
        let response = |object: &str| {
            r#"{"object": OBJECT, "orbit": {"epoch": 2461000.5, "elements": [{"name": "a", "value": 2.0}, {"name": "e", "value": 0.1}, {"name": "i", "value": 5.0}, {"name": "n", "value": 0.3}]}}"#
                .replace("OBJECT", object)
        };
        for (object, requested, expected) in [
            (r#"{"spkid":"20000001"}"#, 1, Some(2_000_001)),
            (r#"{"spkid":2000001}"#, 1, Some(2_000_001)),
            (r#"{"spkid":"2000001"}"#, 0, Some(2_000_001)),
            (r#"{"spkid":"20000373"}"#, 373, Some(2_000_373)),
            (r#"{"spkid":2000373}"#, 373, Some(2_000_373)),
            (r#"{"fullname":"Ceres"}"#, 1, Some(2_000_001)),
            (r#"{"fullname":"Ceres"}"#, 0, None),
            (r#"{"spkid":"20000002"}"#, 1, None),
            (r#"{"spkid":"not-an-identity"}"#, 1, None),
            (r#"{"spkid":"20000001.5"}"#, 1, None),
            (r#"{"spkid":null}"#, 1, None),
            (r#"{"spkid":"20000000"}"#, 1, None),
            (r#"{}"#, i64::MAX, None),
        ] {
            let elements = parse_sbdb_elements_response(&response(object), requested)
                .expect("identity metadata must not alter numerical parsing");
            assert_eq!(
                elements.physical_spice_id(),
                expected,
                "object={object}, request={requested}"
            );
            assert_eq!(elements.epoch().to_bits(), 2461000.5_f64.to_bits());
            assert_eq!(elements.a().to_bits(), 2.0_f64.to_bits());
            assert_eq!(elements.e().to_bits(), 0.1_f64.to_bits());
            assert_eq!(elements.i().to_bits(), 5.0_f64.to_bits());
            assert_eq!(elements.n().to_bits(), 0.3_f64.to_bits());
            assert_eq!(
                (elements.omega(), elements.node(), elements.m0()),
                (0.0, 0.0, 0.0)
            );
            let scalar_copy = OrbitalElements::new(
                elements.name().into(),
                elements.epoch(),
                elements.a(),
                elements.e(),
                elements.i(),
                elements.omega(),
                elements.node(),
                elements.m0(),
                elements.n(),
            )
            .unwrap();
            assert_eq!(scalar_copy.physical_spice_id(), None);
        }
    }

    #[test]
    fn sbdb_name_parser_shapes() {
        assert_eq!(
            parse_sbdb_name_response(r#"{"object": {"spkid": "20099942"}}"#),
            Some(99942)
        );
        assert_eq!(
            parse_sbdb_name_response(r#"{"object": {"spkid": 2000001}}"#),
            Some(1)
        );
        assert_eq!(
            parse_sbdb_name_response(r#"{"object": {"spkid": "2002060"}}"#),
            Some(2060)
        );
        assert_eq!(
            parse_sbdb_name_response(r#"{"object": {"des": "433"}}"#),
            Some(433)
        );
        assert_eq!(
            parse_sbdb_name_response(r#"{"object": {"fullname": "1 Ceres"}}"#),
            Some(1)
        );
        assert!(parse_sbdb_name_response(r#"{"error": "nope"}"#).is_none());
        assert!(parse_sbdb_name_response(r#"{"object": {"spkid": "20000000"}}"#).is_none());
        assert!(parse_sbdb_name_response(r#"{"object": {}}"#).is_none());
    }

    #[test]
    fn spkid_decode_vectors() {
        assert_eq!(decode_sbdb_spkid("20000001"), Some(1));
        assert_eq!(decode_sbdb_spkid("2002060"), Some(2060));
        assert_eq!(decode_sbdb_spkid("20000000"), None);
        assert_eq!(decode_sbdb_spkid("0"), None);
        assert_eq!(decode_sbdb_spkid("abc"), None);
    }

    #[test]
    fn sbdb_fetch_cache_and_sealed_legs() {
        use crate::context::EffectiveNetworkPolicy;
        let mut cache = SbdbElementsCache::new();
        cache.insert(433, sample_elements());
        // Cache hit answers without socket work, even sealed.
        let hit = fetch_sbdb_elements(433, 30.0, EffectiveNetworkPolicy::Sealed, &mut cache)
            .expect("cache hit");
        assert!(hit.is_some());
        // Non-positive numbers never touch the network.
        let mut fresh = SbdbElementsCache::new();
        assert!(
            fetch_sbdb_elements(0, 30.0, EffectiveNetworkPolicy::Sealed, &mut fresh)
                .expect("no socket")
                .is_none()
        );
        // Sealed uncached numbers refuse before any socket work.
        let err =
            fetch_sbdb_elements(1, 30.0, EffectiveNetworkPolicy::Sealed, &mut fresh).unwrap_err();
        assert!(matches!(err, Error::NetworkSealed { .. }));
        assert!(fresh.is_empty());
        // Bad timeouts are argument errors, not silent defaults.
        assert!(
            fetch_sbdb_elements(1, f64::NAN, EffectiveNetworkPolicy::Allow, &mut fresh).is_err()
        );
    }

    #[test]
    fn utc_date_labels_and_naif_chain() {
        // J2000 TT renders on the 2000-01-01 UTC date (TT leads UTC by
        // ~64 s there, inside the same civil day).
        assert_eq!(
            tt_jd_to_utc_date(2_451_545.0),
            Some("2000-01-01".to_string())
        );
        assert!(tt_jd_to_utc_date(f64::NAN).is_none());
        // File-first NAIF chain: catalog number, lone small-body,
        // lone non-center target, legacy name deduction, then refusal.
        assert_eq!(
            deduce_downloaded_naif_id(&[2_000_001, 2002060], "2060", Some(2060)),
            Some(2002060)
        );
        assert_eq!(
            deduce_downloaded_naif_id(&[2002060], "2060", None),
            Some(2002060)
        );
        assert_eq!(deduce_downloaded_naif_id(&[499], "Mars", None), Some(499));
        assert_eq!(
            deduce_downloaded_naif_id(&[0, 10, 499, 599], "Mars", None),
            None
        );
        // Legacy name deduction is the last resort: a bare number still yields the
        // legacy id even when the kernel serves nothing.
        assert_eq!(deduce_downloaded_naif_id(&[], "2060", None), Some(2002060));
        assert_eq!(deduce_downloaded_naif_id(&[], "Nope", None), None);
    }

    #[test]
    fn auto_download_refusal_legs_without_network() {
        use crate::context::EffectiveNetworkPolicy;
        use crate::spk::SpkRegistry;
        let dir = std::env::temp_dir().join(format!("libeph-minornet-{}", std::process::id()));
        let mut registry = SpkRegistry::new();
        // Unsupported bodies never touch the network, even allowed.
        assert!(
            auto_download_asteroid_spk_to(
                999_999_999,
                None,
                None,
                false,
                0.0,
                &dir,
                &mut registry,
                EffectiveNetworkPolicy::Allow,
                30.0
            )
            .expect("unsupported")
            .is_none()
        );
        // Entirely out-of-window requests answer None before socket work.
        assert!(
            auto_download_asteroid_spk_to(
                C_CERES,
                Some(1_000_000.0),
                Some(1_000_100.0),
                false,
                0.0,
                &dir,
                &mut registry,
                EffectiveNetworkPolicy::Allow,
                30.0
            )
            .expect("out of window")
            .is_none()
        );
        // Bad padding is an argument error, not a silent clamp.
        assert!(
            auto_download_asteroid_spk_to(
                C_CERES,
                None,
                None,
                false,
                -1.0,
                &dir,
                &mut registry,
                EffectiveNetworkPolicy::Allow,
                30.0
            )
            .is_err()
        );
        // Ensure on an unregistered body under a sealed policy refuses
        // (never a silent `false` from a swallowed refusal).
        let err = ensure_major_asteroid_spk_to(
            C_CERES,
            None,
            &dir,
            &mut registry,
            EffectiveNetworkPolicy::Sealed,
            30.0,
        )
        .unwrap_err();
        assert!(matches!(err, Error::NetworkSealed { .. }));
    }

    #[test]
    fn auto_download_membership_short_circuit() {
        use crate::context::EffectiveNetworkPolicy;
        use crate::spk::SpkRegistry;
        let dir = std::env::temp_dir().join(format!("libeph-minorreg-{}", std::process::id()));
        let mut registry = SpkRegistry::new();
        // The mapping below only exercises the membership short-circuit
        // (the stored path answers before any ephemeris evaluation, so
        // the target is never propagated).
        registry
            .register_image("k.bsp", crate::spk::moon_test_kernel_bytes())
            .expect("image");
        registry
            .register_body(C_CERES, "k.bsp", 3)
            .expect("earth-segment mapping");
        let got = auto_download_asteroid_spk_to(
            C_CERES,
            None,
            None,
            false,
            0.0,
            &dir,
            &mut registry,
            EffectiveNetworkPolicy::Sealed,
            30.0,
        )
        .expect("membership answers without socket work");
        assert_eq!(got, Some(std::path::PathBuf::from("k.bsp")));
        assert!(
            ensure_major_asteroid_spk_to(
                C_CERES,
                None,
                &dir,
                &mut registry,
                EffectiveNetworkPolicy::Sealed,
                30.0,
            )
            .expect("ensure on registered body")
        );
    }

    #[test]
    fn asteroid_by_number_unknown_body_without_provider() {
        // An SBDB-only number with an empty provider is UnknownBody
        // ; sealed refusals still propagate.
        let empty = |_: i64| -> Result<Option<OrbitalElements>> { Ok(None) };
        let err = calc_asteroid_by_number(999_999_937, 2_461_000.5, false, None, Some(&empty), 0.0)
            .unwrap_err();
        assert!(matches!(err, Error::UnknownBody { .. }));
        let sealed = |_: i64| -> Result<Option<OrbitalElements>> {
            Err(Error::NetworkSealed {
                purpose: "probe".to_string(),
            })
        };
        let err =
            calc_asteroid_by_number(999_999_937, 2_461_000.5, false, None, Some(&sealed), 0.0)
                .unwrap_err();
        assert!(matches!(err, Error::NetworkSealed { .. }));
        // Local numbers serve without any provider.
        let (lon, lat, dist) =
            calc_asteroid_by_number(1, 2_451_545.0, true, None, None, 0.0).expect("ceres");
        assert!((0.0..360.0).contains(&lon));
        assert!((-90.0..=90.0).contains(&lat));
        assert!(dist > 0.0);
    }

    #[test]
    fn caches_and_local_lookup() {
        let mut names = AsteroidNameCache::new();
        assert_eq!(names.get_local("Ceres"), Some(1));
        assert_eq!(names.get_local("  EROS  "), Some(433));
        assert_eq!(names.get_local("europa"), Some(52));
        assert_eq!(names.get_local("No Such Body"), None);
        assert_eq!(names.get_local("   "), None);
        // Seeding is idempotent: a remote insert survives a second seed.
        names.insert("Dinky", 999001);
        names.seed_local();
        assert_eq!(names.get_local("dinky"), Some(999001));
        let cleared = names.clear();
        assert_eq!(cleared, LOCAL_ASTEROID_NAMES.len() + 1);
        assert!(names.is_empty());

        let mut elements = SbdbElementsCache::new();
        assert!(elements.get(433).is_none());
        elements.insert(433, sample_elements());
        assert!(elements.get(433).is_some());
        assert_eq!(elements.clear(), 1);

        // Remote insert is cached under the normalized name.
        let mut cache = AsteroidNameCache::new();
        let out =
            get_asteroid_number("Dinky", &mut cache, Some(&|_| Ok(Some(999001)))).expect("lookup");
        assert_eq!(out, Some(999001));
        let again = get_asteroid_number("dinky", &mut cache, None).expect("cached");
        assert_eq!(again, Some(999001));
        // Transport failures (other than sealed) become misses.
        let miss = get_asteroid_number(
            "X",
            &mut cache,
            Some(&|_| {
                Err(Error::HttpTransport {
                    purpose: "t".to_string(),
                    detail: "d".to_string(),
                })
            }),
        )
        .expect("swallowed");
        assert_eq!(miss, None);
        assert!(
            get_asteroid_number("X", &mut cache, None)
                .expect("none")
                .is_none()
        );
    }

    #[test]
    fn sealed_name_lookup_propagates() {
        let mut cache = AsteroidNameCache::new();
        let sealed = Error::NetworkSealed {
            purpose: "SBDB lookup".to_string(),
        };
        let out = get_asteroid_number("Zed", &mut cache, Some(&|_| Err(sealed.clone())));
        assert!(matches!(out, Err(Error::NetworkSealed { .. })));
    }

    #[test]
    fn spk_window_clamp_shapes() {
        assert_eq!(
            clamp_spk_window(2451545.0, 2461000.5, 0.0).expect("inside"),
            Some((2451545.0, 2461000.5))
        );
        // Partially outside clamps to the padded range.
        assert_eq!(
            clamp_spk_window(2000000.0, 2461000.5, 0.0).expect("clamped"),
            Some((HORIZONS_SPK_JD_MIN, 2461000.5))
        );
        // Entirely outside refuses.
        assert_eq!(
            clamp_spk_window(1000000.0, 2000000.0, 0.0).expect("outside"),
            None
        );
        // Padding extends the range.
        assert_eq!(
            clamp_spk_window(HORIZONS_SPK_JD_MIN - 5.0, HORIZONS_SPK_JD_MAX + 5.0, 10.0)
                .expect("padded"),
            Some((HORIZONS_SPK_JD_MIN - 5.0, HORIZONS_SPK_JD_MAX + 5.0))
        );
        assert!(matches!(
            clamp_spk_window(2451545.0, 2461000.5, -1.0),
            Err(Error::InvalidPadding { .. })
        ));
        assert!(matches!(
            clamp_spk_window(2451545.0, 2461000.5, f64::NAN),
            Err(Error::InvalidPadding { .. })
        ));
    }

    #[test]
    fn spk_targets_and_registry() {
        let (number, horizons, naif, name) = spk_download_target(C_CERES).expect("ceres target");
        assert_eq!((number, horizons, naif), (1, "Ceres;", 2000001));
        assert_eq!(name, "Ceres");
        assert!(spk_download_target(999999).is_none());
        let majors = list_major_asteroids();
        assert_eq!(
            majors.iter().map(|(id, _)| *id).collect::<Vec<_>>(),
            vec![C_CERES, C_PALLAS, C_JUNO, C_VESTA, C_CHIRON]
        );
        assert!(is_spk_available_for_body(&|_| true, C_CERES));
        assert!(!is_spk_available_for_body(&|_| false, C_CERES));
        // Reused SPK-07 gates agree with the model rows.
        assert!(crate::spk::is_major_asteroid(C_CERES));
        assert!(!crate::spk::is_major_asteroid(C_ERIS));
        assert!(!crate::spk::is_in_jpl_major_body_index(C_CERES));
    }

    #[test]
    fn heliocentric_entries_and_fallbacks() {
        // Unknown body refuses on both entries.
        assert!(matches!(
            calc_minor_body_heliocentric(999999, 2461000.5, false, None, 0.0),
            Err(Error::UnknownBody { .. })
        ));
        assert!(matches!(
            calc_asteroid_by_number(999999, 2461000.5, false, None, None, 0.0),
            Err(Error::UnknownBody { .. })
        ));
        // Bad catalog numbers refuse before any provider work.
        assert!(matches!(
            calc_asteroid_by_number(0, 2461000.5, false, None, None, 0.0),
            Err(Error::InvalidBody { .. })
        ));
        assert!(matches!(
            resolve_asteroid_number(-3),
            Err(Error::InvalidBody { .. })
        ));
        // Table + special routing agree with the direct entry.
        let direct = calc_minor_body_heliocentric(C_CERES, 2461000.5, false, None, 0.0)
            .expect("ceres direct");
        let via_number =
            calc_asteroid_by_number(1, 2461000.5, false, None, None, 0.0).expect("ceres map");
        assert_eq!(direct, via_number);
        // SPK hit short-circuits.
        let spk_hit = calc_minor_body_heliocentric(
            C_CERES,
            2461000.5,
            true,
            Some(&|_, _| Ok(Some((10.0, 20.0, 3.0)))),
            0.0,
        )
        .expect("spk hit");
        assert_eq!(spk_hit, (10.0, 20.0, 3.0));
        // Range miss falls through to the Keplerian path.
        let fell_through = calc_minor_body_heliocentric(
            C_CERES,
            2461000.5,
            true,
            Some(&|_, _| {
                Err(Error::EphemerisRange {
                    requested_jd: 0.0,
                    start_jd: 0.0,
                    end_jd: 1.0,
                    body: None,
                })
            }),
            0.0,
        )
        .expect("fallthrough");
        assert_eq!(fell_through, direct);
        // Fatal provider failures propagate.
        assert!(matches!(
            calc_minor_body_heliocentric(
                C_CERES,
                2461000.5,
                true,
                Some(&|_, _| Err(Error::SpkRequired { body: C_CERES })),
                0.0,
            ),
            Err(Error::SpkRequired { .. })
        ));
        // SBDB leg composes like the direct path with body_id = None.
        let elements = minor_body_elements(C_VESTA).expect("vesta");
        let via_sbdb = calc_asteroid_by_number(
            500000,
            2461000.5,
            false,
            None,
            Some(&|_| Ok(Some(elements.clone()))),
            0.0,
        )
        .expect("sbdb leg");
        let [x, y, z] = calc_minor_body_position(&elements, 2461000.5, true, None).expect("raw");
        let (lon, lat, r) = cartesian_to_spherical(x, y, z).expect("shape");
        let (lon_p, lat_p) =
            crate::hypothetical::precess_ecliptic(lon, lat, J2000_EPOCH, 2461000.5)
                .expect("precess");
        assert_eq!(via_sbdb, (lon_p, lat_p, r));
        // Sealed SBDB refusal propagates.
        assert!(matches!(
            calc_asteroid_by_number(
                500000,
                2461000.5,
                false,
                None,
                Some(&|_| Err(Error::NetworkSealed {
                    purpose: "x".to_string()
                })),
                0.0,
            ),
            Err(Error::NetworkSealed { .. })
        ));
        // Non-finite epochs refuse.
        assert!(calc_minor_body_heliocentric(C_CERES, f64::NAN, false, None, 0.0).is_err());
        assert!(calc_minor_body_position(&elements, f64::INFINITY, true, None).is_err());
    }

    #[test]
    fn resolve_routing_shapes() {
        // Eros (433) is predefined: table wins.
        assert_eq!(
            resolve_asteroid_number(433).expect("eros"),
            Some(433 + AST_OFFSET)
        );
        // Ceres number hits the predefined table (special map is the fallback).
        assert_eq!(resolve_asteroid_number(1).expect("ceres"), Some(C_CERES));
        assert_eq!(
            resolve_asteroid_number(2060).expect("chiron"),
            Some(C_CHIRON)
        );
        // Unlisted, non-special numbers leave only the SBDB leg.
        assert_eq!(resolve_asteroid_number(1234567).expect("unlisted"), None);
    }

    #[test]
    fn forced_vectors_match_release_coupling_order() {
        // Forced secular vectors `(g, h, k, s, p, q)` with `g`/`s` in
        // rad/day and raw J2000 propagation `(x, y, z)` in AU for Ceres at
        // one second before JD 2460400.0 TT (`2460400.0 - 1/86400`). The
        // planet-vector coupling multiplies the formed planet vectors
        // (`nu_j * h_j`, `sigma_j * p_j`); the fused triple product rounds
        // a planet term 1 ulp off and the 1-second nod_aps stencil
        // amplifies the slop ~43000 x on this near-circular orbit.
        let jd = 2460400.0 - 1.0 / 86400.0;
        let (primary, _, _) = epoch_elements_blend(C_CERES, jd).expect("blend");
        let (g, h, k, s, p, q) = forced_elements(&primary, jd);
        assert_eq!(g.to_bits(), 0x3ec96b20a4fd7e1a);
        assert_eq!(h.to_bits(), 0x3f7154ba15a478c9);
        assert_eq!(k.to_bits(), 0x3f8c89373f0d975f);
        assert_eq!(s.to_bits(), 0x3eb8c9ea6e8a7e77);
        assert_eq!(p.to_bits(), 0x3f8826327e0c659e);
        assert_eq!(q.to_bits(), 0xbf6421786bc817b0);
        let xyz = calc_minor_body_position(&primary, jd, true, Some(C_CERES)).expect("raw");
        assert_eq!(xyz[0].to_bits(), 0xbfcfcb91f88e6a19);
        assert_eq!(xyz[1].to_bits(), 0xc0068ec469016fef);
        assert_eq!(xyz[2].to_bits(), 0xbfa5eebc8c0564f2);
    }
}

#[cfg(test)]
mod nbody_source_record_tests {
    use super::*;
    #[test]
    fn source_record_retains_ieee_without_changing_checked_constructor() {
        let nan = f64::from_bits(0x7ff8_0000_0000_0123);
        let row = OrbitalElements::from_nbody_source_record(
            "\0é".into(),
            2451545.0,
            f64::NEG_INFINITY,
            nan,
            -0.0,
            40.0,
            70.0,
            11.0,
            0.25,
        );
        assert_eq!(row.name(), "\0é");
        assert_eq!(row.a().to_bits(), f64::NEG_INFINITY.to_bits());
        assert_eq!(row.e().to_bits(), nan.to_bits());
        assert_eq!(row.i().to_bits(), (-0.0_f64).to_bits());
        assert_eq!(row.n().to_bits(), 0.25_f64.to_bits());
        assert!(
            OrbitalElements::new(
                "checked".into(),
                2451545.0,
                f64::NEG_INFINITY,
                0.1,
                12.,
                40.,
                70.,
                11.,
                0.25
            )
            .is_err()
        );
    }
}
