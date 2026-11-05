//! Planetary nodes, apsides, osculating elements and distance extrema.
//!
//! Mean elements use Meeus Table 31.A after Simon et al. (1994). Osculating
//! quantities reduce explicit Cartesian states with the selected central mass;
//! Earth-system geometry composes the Earth–Moon barycenter. Node and apsis
//! positions follow the same output frame, observer and sidereal selections as
//! their declared entry. Lunar points combine mean theory, stored channels and
//! same-epoch two-body geometry.
//!
//! Epochs are TT Julian days unless an entry is explicitly UT; UT conversion
//! uses its bound clock. Distances are AU, rates AU/day or degrees/day, and
//! spherical output angles degrees. Low-level conic angles use radians where
//! specified by their entry points. State tuples contain longitude, latitude, range and
//! their three rates. Native geometric derivatives retain actual polynomial
//! ownership and the full epoch precision.
//!
//! A calculation requires the selected provider, coverage, observer and any
//! sidereal anchor. Missing inputs, unsupported objects, invalid domains and
//! nonfinite checked arguments produce typed errors. Elliptic extrema and
//! osculating reductions retain their stated two-body/model limitations; they
//! are not replacements for an N-body trajectory.
//!
//! # Example: mean elements and a point on their reference conic
//!
//! ```
//! use libephemeris_core::{constants::EARTH, orbital};
//!
//! let [a, e, inclination, node, perihelion] =
//!     orbital::mean_orbital_elements(EARTH, 2451545.0).expect("Earth mean elements");
//! let shape = orbital::mean_shape_from_elements(a, e, inclination, node, perihelion);
//! // True anomaly is radians; the vector is AU in the mean ecliptic of date.
//! let (position_au, radius_au) = orbital::orbit_point_3d(shape, 0.0);
//! println!("Perihelion: {position_au:?}, radius {radius_au} AU");
//! ```

#[path = "orbital/elements.rs"]
mod elements;
#[path = "orbital/extrema.rs"]
mod extrema;
#[path = "orbital/geometry.rs"]
mod geometry;
#[path = "orbital/projection.rs"]
mod projection;

// Keep the consumed facade paths stable while each private child owns a
// numerical responsibility. The coefficient helpers remain crate-only for
// the perturber models; they are not part of the public orbital API.
pub use elements::{
    ecliptic_spherical_to_cartesian_state, mean_orbital_elements,
    orbital_elements_from_ecliptic_state, osculating_semi_major_axis,
};
pub(crate) use elements::{mean_element_coeffs, poly_asc};
pub use extrema::{
    earth_extremes_from_emb_elements, ellipse_position_au, geocentric_two_ellipse_extremes,
    helio_distance_extremes, orbit_ellipse_from_elements, sun_extremes_from_emb_elements,
    two_ellipse_min_max,
};
use geometry::orbit_node_3d;
pub use geometry::{
    focal_point_3d, mean_shape_from_elements, moon_nodaps_points, orbit_point_3d,
    osculating_shape_from_ecliptic,
};
pub use projection::{
    apply_nodaps_j2000_rule, nodaps_dir_to_j2000, nodaps_point_speeds,
    nodaps_ray_crosses_solar_disk, nodaps_vector_to_j2000,
};

use crate::constants;
use crate::errors::{Error, Result};
use crate::positions::{
    C_LIGHT_AU_DAY, PositionResult, SealedServing, SiteEllipsoid, StateProvider, apply_aberration,
    apply_output_flags, icrs_to_ecl_j2000_matrix, mat3_vec3, serving_channel_gate,
    serving_compute_tt, topocentric_offset, vondrak_mean_obliquity_deg,
};
use crate::sidereal::{SealedAyanamshaServing, is_fixed_epoch_request, strip_to_base_mode};

// ---------------------------------------------------------------------------
// Central masses, tables and method contracts
// ---------------------------------------------------------------------------

/// Gaussian gravitational constant `k` of the IAU (1976) system
/// : the IAU (1976) System of
/// Astronomical Constants, Explanatory Supplement to the Astronomical
/// Almanac, 3 rd ed., Tables 8.1/15.1.
pub const GAUSS_K: f64 = 0.01720209895;

/// Solar gravitational parameter in AU^3/day^2 (`k^2`; the squaring is one
/// correctly-rounded multiply on both sides, so the spelling is exact).
pub const GM_SUN: f64 = GAUSS_K * GAUSS_K;

/// TT Julian Day of J2000.0.
pub const J2000_JD: f64 = 2451545.0;

/// Reciprocal system-mass ratios Sun/planet-system.
fn sun_mass_ratio(body_id: i32) -> Option<f64> {
    match body_id {
        constants::MERCURY => Some(6023600.0),
        constants::VENUS => Some(408523.71),
        constants::EARTH => Some(328900.5596),
        constants::MARS => Some(3098708.0),
        constants::JUPITER => Some(1047.3486),
        constants::SATURN => Some(3497.898),
        constants::URANUS => Some(22902.98),
        constants::NEPTUNE => Some(19412.24),
        constants::PLUTO => Some(135200000.0),
        _ => None,
    }
}

/// Default two-body parameter `G*(M_sun + M_system)` for a major planet,
/// or `None` where the model table has no entry (the caller then keeps
/// `GM_SUN`, exactly like `mass_ratio` fallback). Ratios are the
/// DE440 planetary masses (Park, Folkner, Williams & Boggs 2021, AJ 161,
/// 105; cf. Astronomical Almanac Table K 7), Earth counting as the
/// Earth+Moon system. Keeping the planet mass matters: for Jupiter it moves
/// the perihelion direction at the ~0.1 deg level.
pub fn planet_system_gm(body_id: i32) -> Option<f64> {
    sun_mass_ratio(body_id).map(|ratio| GM_SUN * (1.0 + 1.0 / ratio))
}

/// Sun-to-Earth (alone) mass ratio, IAU 2009 nominal: Luzum et al. 2011, Celestial Mechanics and
/// Dynamical Astronomy 110, 293. For geocentric osculating reductions only,
/// never for heliocentric planet orbits.
pub const SUN_EARTH_MASS_RATIO: f64 = 332946.0487;

/// Earth+Moon two-body parameter for the Moon's geocentric osculating
/// orbit: the
/// about-Earth two-body mu. The Earth-alone ratio would bias the
/// semi-major axis by about +1.2% here.
pub fn moon_system_gm() -> f64 {
    GM_SUN / 328900.5596
}

/// Earth-only parameter for the Moon's derived period slots. It leaves a ~8 e-4 relative offset on the period family (down
/// from 5.3 e-3 with GM(Earth+Moon)); no externally realized geocentric GM
/// is published for comparison.
pub fn moon_central_gm() -> f64 {
    GM_SUN / SUN_EARTH_MASS_RATIO
}

/// Semi-major-axis classifier for barycentric re-referencing: only orbits beyond Jupiter respond to
/// `FLG_BARYCTR`. Any threshold in (5.21, 9.53) classifies identically;
/// 6.0 is the model project choice.
pub const TRANS_JOVIAN_A_AU: f64 = 6.0;

/// IAU 2015 nominal solar radius in AU: the finite-disk gate for the PPN
/// point-mass deflection.
pub const NOMINAL_SOLAR_RADIUS_AU: f64 = 6.957e8 / 149_597_870_700.0;

/// Spacing of the fourth-order node/apse derivative stencil, in TT days.
///
/// The dyadic value is exactly representable when added to supported Julian
/// epochs. Its 22.5-minute spacing limits roundoff amplification in nearly
/// stationary angular coordinates; the fourth-order stencil suppresses the
/// truncation error that a wider two-point difference would introduce.
pub const NODAPS_SPEED_STEP_DAYS: f64 = 1.0 / 64.0;

/// Textbook mean lunar distance in AU: the mean-ellipse identity
/// `r_peri = 2a - r_apog`.
pub const MEAN_LUNAR_DISTANCE_AU: f64 = 0.002569555;

/// Major-planet order for the Astronomical-Almanac Jacobi hierarchy
const ORBEL_AA_ORDER: [i32; 9] = [
    constants::MERCURY,
    constants::VENUS,
    constants::EARTH,
    constants::MARS,
    constants::JUPITER,
    constants::SATURN,
    constants::URANUS,
    constants::NEPTUNE,
    constants::PLUTO,
];

/// Nominal semi-major axes (AU) ordering the hierarchy: only the ordering matters, not precision.
fn orbel_aa_nominal_a(body_id: i32) -> Option<f64> {
    match body_id {
        constants::MERCURY => Some(0.387),
        constants::VENUS => Some(0.723),
        constants::EARTH => Some(1.000),
        constants::MARS => Some(1.524),
        constants::JUPITER => Some(5.204),
        constants::SATURN => Some(9.582),
        constants::URANUS => Some(19.19),
        constants::NEPTUNE => Some(30.07),
        constants::PLUTO => Some(39.48),
        _ => None,
    }
}

/// Summed reciprocal system-mass ratios through the body's own system
/// : the `FLG_ORBEL_AA` central-mass
/// excess for a major planet. Accumulates outward in hierarchy order.
/// Hierarchical Jacobi
/// parameter `GM_eff(k) = G*(M_sun + m_1 + ... + m_k)`: the Sun plus every
/// planetary system at or interior to the body's own orbit, self included;
/// a planet counts with its system mass (planet plus satellites), Earth
/// entering as the Earth+Moon barycentre. Only the osculating geometry
/// responds; the derived period family keeps pure solar GM.
pub fn orbel_aa_cumulative_through_self(planet: i32) -> Option<f64> {
    let mut total = 0.0;
    for pid in ORBEL_AA_ORDER {
        let ratio = sun_mass_ratio(pid)?;
        total += 1.0 / ratio;
        if pid == planet {
            return Some(total);
        }
    }
    None
}

/// Summed reciprocal system-mass ratios strictly interior to `a_au`
/// : the `FLG_ORBEL_AA` excess for
/// a minor body (restricted-problem Jacobi parameter, no self term): the
/// Sun plus every major planet strictly interior to its orbit. The body's
/// own negligible mass adds no term and minor bodies are not hierarchy
/// members. The nominal axes only order the nine systems; their precision
/// is immaterial.
pub fn orbel_aa_interior_mass_sum(a_au: f64) -> f64 {
    let mut total = 0.0;
    for pid in ORBEL_AA_ORDER {
        if let Some(nominal) = orbel_aa_nominal_a(pid) {
            if nominal < a_au {
                if let Some(ratio) = sun_mass_ratio(pid) {
                    total += 1.0 / ratio;
                }
            }
        }
    }
    total
}

/// Mean-vs-osculating model selection from the `method` low byte
/// : mean when the low byte is
/// zero or `NODBIT_MEAN` is set, osculating otherwise. Bit 8 (`256`,
/// `NODBIT_FOPOINT`) never influences the choice. Accepts every integer;
/// Bitwise `&` on `i32` performs low-byte masking
/// (two's complement).
pub fn prefer_mean(method: i32) -> bool {
    (method & 0xFF) == 0 || (method & constants::NODBIT_MEAN) != 0
}

/// Bodies with no orbital elements for the elements entries: lunar nodes/apogees plus the
/// interpolated apsides. The Sun is handled separately by the caller.
fn is_orbital_invalid_object(body_id: i32) -> bool {
    matches!(
        body_id,
        constants::MEAN_NODE
            | constants::TRUE_NODE
            | constants::MEAN_APOG
            | constants::OSCU_APOG
            | constants::INTP_APOG
            | constants::INTP_PERG
    )
}

/// Reject ids with no orbital elements: `object N not valid` for negatives, the
/// lunar nodes/apogees, the interpolated apsides (and the Sun unless
/// `sun_valid`), `illegal planet number N` for 23-39, and `object N not
/// valid` for the sourceless block between the last fictitious body and
/// the asteroid offset. That block has no element source; without the
/// guard the calculator would leak a zero-initialized tuple reading as a
/// real 0-AU orbit. `entry` names the public caller and reproduces
/// the model message as `"{entry}: {detail}"`.
pub fn validate_orbital_object(entry: &'static str, body_id: i32, sun_valid: bool) -> Result<()> {
    if body_id < 0
        || is_orbital_invalid_object(body_id)
        || (body_id == constants::SUN && !sun_valid)
    {
        return Err(Error::InvalidOrbitalObject {
            entry,
            detail: format!("error in get_orbital_elements(): object {body_id} not valid"),
        });
    }
    if (23..=39).contains(&body_id) {
        return Err(Error::InvalidOrbitalObject {
            entry,
            detail: format!("illegal planet number {body_id}."),
        });
    }
    if constants::WALDEMATH < body_id && body_id < constants::AST_OFFSET {
        return Err(Error::InvalidOrbitalObject {
            entry,
            detail: format!("error in get_orbital_elements(): object {body_id} not valid"),
        });
    }
    Ok(())
}

/// membership for the node/apse entries:
/// Sun, Moon, Mercury-Neptune, Pluto and Earth.
pub fn is_nodaps_planet(body_id: i32) -> bool {
    matches!(
        body_id,
        constants::SUN
            | constants::MOON
            | constants::MERCURY
            | constants::VENUS
            | constants::EARTH
            | constants::MARS
            | constants::JUPITER
            | constants::SATURN
            | constants::URANUS
            | constants::NEPTUNE
            | constants::PLUTO
    )
}

/// membership plus the numbered-asteroid
/// range (`AST_OFFSET < ipl < FIXSTAR_OFFSET`, exclusive on both ends for
/// the node/apse routing).
pub fn is_minor_nodaps_body(body_id: i32) -> bool {
    matches!(
        body_id,
        constants::CHIRON
            | constants::PHOLUS
            | constants::CERES
            | constants::PALLAS
            | constants::JUNO
            | constants::VESTA
    ) || (constants::AST_OFFSET < body_id && body_id < constants::FIXSTAR_OFFSET)
}

// ---------------------------------------------------------------------------
// Published mean orbital elements (Meeus Table 31.A / Simon et al. 1994)
// ---------------------------------------------------------------------------

/// Coefficient bundle behind [`mean_orbital_elements`].
pub(crate) struct MeanCoeffs {
    pub(crate) a: &'static [f64],
    pub(crate) e: &'static [f64],
    pub(crate) i: &'static [f64],
    pub(crate) omega: Option<&'static [f64]>,
    pub(crate) pi: &'static [f64],
}

// ---------------------------------------------------------------------------
// Node/apse orbit-frame geometry
// ---------------------------------------------------------------------------

/// Osculating (or mean) orbit shape in the working ecliptic frame:
/// semi-major axis, eccentricity, inclination, ascending node and
/// argument of perihelion (radians) plus the semi-latus rectum.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OrbitShape {
    a: f64,
    e: f64,
    incl: f64,
    node: f64,
    arg_peri: f64,
    p: f64,
}

// ---------------------------------------------------------------------------
// Working-frame rotation (ICRS <-> ecliptic of date)
// ---------------------------------------------------------------------------

/// Rigid-rotation arm for a projected `FLG_SIDEREAL` node/apse request
/// : the ecliptic
/// set rotates into the target frame, or the equatorial set into the
/// mean equator of t 0.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum NodapsProjectionRoute {
    Ecliptic,
    EquatorialT0,
}

/// Native ICRS-to-ecliptic rotation for one node/apse computation.
///
/// Uses the same direct TT Vondrák pole basis as the main position pipeline.
/// Nutation rotates longitude about the ecliptic pole; obliquity nutation
/// cannot change ecliptic latitude. A single composed basis avoids recovering
/// tiny ecliptic heights by subtracting large intermediate equatorial terms.
#[derive(Debug, Clone, Copy)]
pub struct NodapsFrame {
    matrix: [[f64; 3]; 3],
}

// ---------------------------------------------------------------------------
// Geocentric projection of one orbit-frame point
// ---------------------------------------------------------------------------

/// Observer/frame context projecting one orbit-frame point to output
/// longitude/latitude/distance: the
/// heliocentric point goes geocentric (minus the observer) unless the
/// request is centered; the mirrored Sun negates instead; aberration
/// and the PPN deflection apply under the model flag gates. Fields
/// are private under the construction invariant; build with [`NodapsProjection::new`].
pub struct NodapsProjection<'d> {
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
    optional_deflect_at: Option<&'d NodapsDeflectorAt<'d>>,
}

/// Fourth-order finite-difference support selected from the containing record.
/// The nominal dyadic step is reduced only to remain within the common domain.
fn nodaps_derivative_stencil(epoch: f64, start: f64, end: f64) -> Result<([f64; 4], f64)> {
    let mut step = NODAPS_SPEED_STEP_DAYS;
    loop {
        if epoch - 2.0 * step >= start && epoch + 2.0 * step <= end {
            return Ok(([-2.0, -1.0, 1.0, 2.0], step));
        }
        if epoch + 4.0 * step <= end && epoch >= start {
            return Ok(([1.0, 2.0, 3.0, 4.0], step));
        }
        if epoch - 4.0 * step >= start && epoch <= end {
            return Ok(([-1.0, -2.0, -3.0, -4.0], step));
        }
        step *= 0.5;
        if step == 0.0 || epoch + step == epoch || epoch - step == epoch {
            return Err(Error::InvalidOrbitalObject {
                entry: "nod_aps",
                detail: "no representable derivative stencil within the common polynomial support"
                    .to_string(),
            });
        }
    }
}

/// Differentiate displacements from the central position, avoiding constant
/// coordinate cancellation. Longitude displacements are unwrapped by the caller.
fn nodaps_stencil_derivative(delta: [f64; 4], offsets: [f64; 4], step: f64) -> f64 {
    let weights = if offsets[0] == -2.0 {
        [1.0, -8.0, 8.0, -1.0]
    } else if offsets[0] > 0.0 {
        [48.0, -36.0, 16.0, -3.0]
    } else {
        [-48.0, 36.0, -16.0, 3.0]
    };
    let mut sum = 0.0;
    let mut correction = 0.0;
    for (value, weight) in delta.into_iter().zip(weights) {
        let term = value * weight;
        let next = sum + term;
        correction += if sum.abs() >= term.abs() {
            (sum - next) + term
        } else {
            (term - next) + sum
        };
        sum = next;
    }
    (sum + correction) / (12.0 * step)
}

/// Undefined-point sentinel for the interpolated lunar apsides
/// : NaN position slots everywhere; the speed
/// slots stay 0.0 on a plain request but propagate NaN under an
/// explicit `FLG_SPEED` request.
pub fn intp_undefined_points(speed: bool) -> [PositionResult; 4] {
    let spd = if speed { f64::NAN } else { 0.0 };
    let pt = [f64::NAN, f64::NAN, f64::NAN, spd, spd, spd];
    [pt, pt, pt, pt]
}

// ---------------------------------------------------------------------------
// Sealed serving
// ---------------------------------------------------------------------------

/// Apply the output-flag shaping layer to all four points (with the topo arm running at the entry
/// before this layer (the calculation order: topo, then equatorial, then
/// reshaping last).
fn shape_nodaps_output(pts: [PositionResult; 4], flags: i32) -> [PositionResult; 4] {
    [
        apply_output_flags(pts[0], flags),
        apply_output_flags(pts[1], flags),
        apply_output_flags(pts[2], flags),
        apply_output_flags(pts[3], flags),
    ]
}

/// First-order diurnal Bradley aberration for the native node/apse site tail.
///
/// Position is ICRS/output-frame AU and relative velocity is AU/day. Normalize
/// the direction twice as in the source Cartesian helper, add beta minus its
/// radial projection, renormalize and restore distance. Magnitudes below
/// 1 e-30 preserve the input. Undefined points propagate NaNs through this
/// private reducer; they must still acquire their finite observer support.
/// This helper adds no provider reads or relativistic light-time branch.
// Source: astrometry.apply_aberration_to_position/_annual_aberration_cartesian.
fn nodaps_diurnal_aberration(position: [f64; 3], velocity: [f64; 3]) -> [f64; 3] {
    let [x, y, z] = position;
    let distance = (x * x + y * y + z * z).sqrt();
    if distance < 1e-30 {
        return position;
    }
    let [ux, uy, uz] = [x / distance, y / distance, z / distance];
    let direction_norm = (ux * ux + uy * uy + uz * uz).sqrt();
    let displacement = if direction_norm < 1e-30 {
        [0.0; 3]
    } else {
        let [rx, ry, rz] = [
            ux / direction_norm,
            uy / direction_norm,
            uz / direction_norm,
        ];
        let [vx, vy, vz] = velocity.map(|v| v / C_LIGHT_AU_DAY);
        let dot = rx * vx + ry * vy + rz * vz;
        [vx - dot * rx, vy - dot * ry, vz - dot * rz]
    };
    let [ax, ay, az] = [
        ux + displacement[0],
        uy + displacement[1],
        uz + displacement[2],
    ];
    let norm = (ax * ax + ay * ay + az * az).sqrt();
    if norm < 1e-30 {
        return position;
    }
    [
        ax / norm * distance,
        ay / norm * distance,
        az / norm * distance,
    ]
}

/// Output selection for the shared conic tail. Constructed from the request
/// flags before acquisition so callers cannot supply conflicting flag views.
struct NodapsPointRequest {
    body_id: i32,
    jd_tt: f64,
    method: i32,
    flags: i32,
    helio: bool,
    bary: bool,
    want_sidereal: bool,
    want_equatorial: bool,
}

impl NodapsPointRequest {
    fn new(body_id: i32, jd_tt: f64, method: i32, flags: i32) -> Self {
        let helio = (flags & constants::FLG_HELCTR) != 0;
        // Node/apse center precedence differs from calc: heliocentric wins
        // when both center flags are set.
        let bary = (flags & constants::FLG_BARYCTR) != 0 && !helio;
        let want_sidereal = (flags & constants::FLG_SIDEREAL) != 0;
        let want_equatorial = (flags & constants::FLG_EQUATORIAL) != 0;
        Self {
            body_id,
            jd_tt,
            method,
            flags,
            helio,
            bary,
            want_sidereal,
            want_equatorial,
        }
    }
}

/// Conic fit and its origin/frame interpretation. Mean fits remain
/// heliocentric; osculating fits carry the caller's resolved barycentric mode.
struct NodapsOrbitFit {
    shape: OrbitShape,
    mean_frame: bool,
    mirror_sun: bool,
    bar_mode: bool,
}

/// Already acquired support for conic projection. These are stack-owned
/// vectors, not providers or retained callbacks. Each caller assembles this
/// record only after its original support reads and frame rotations.
struct NodapsAcquiredObserver {
    frame: NodapsFrame,
    nut_eff: Option<(f64, f64)>,
    sun_p: [f64; 3],
    earth_p: [f64; 3],
    earth_helio_ecl: [f64; 3],
    earth_bary_ecl: [f64; 3],
    earth_vel_ecl: [f64; 3],
}

#[derive(Clone, Copy)]
struct NodapsSmoothSupport {
    epoch: f64,
    start: f64,
    end: f64,
}

/// Shared node/apse and orbital-element reduction over borrowed geometric
/// ICRS states and explicit epoch/acquisition callbacks. The stored route uses
/// caller Delta-T; optional Native geometry retains original split Time and
/// independent observation clocks. Fields are private; no globals or I/O.
/// Minor/fictitious reductions can borrow a public TT calculation feed through
/// [`Self::with_calc_tt`]. The independent observation timescale for public
/// orbital UT is selected separately through [`Self::with_orbital_ut1_to_tt`].
///
/// Nutation for the of-date node/apse working frame is the native IAU
/// 2006/2000 A pair in radians, evaluated at each central and stencil TT
/// epoch; see the module docs. `FLG_TOPOCTR` without an explicit site
/// is refused when parallax applies; centered views suppress it. The
/// site-carrying
/// [`SealedOrbitalServing::nod_aps_topo_tt`] entry serves it.
/// `FLG_SIDEREAL` needs the
/// attached ayanamsha serving ([`SealedOrbitalServing::with_ayanamsha`])
/// and runs the plain mean-ayanamsha arm, or the rigid frame projection
/// for fixed-epoch / `SIDBIT_ECL_T0` / `SIDBIT_SSY_PLANE` selections;
/// otherwise it keeps the same refusal, which is also the production
/// binding behavior.
pub struct SealedOrbitalServing<'a, P> {
    serving: &'a SealedServing<'a, P>,
    deltat_at: &'a dyn Fn(f64) -> Result<f64>,
    ayanamsha: Option<SealedAyanamshaServing<'a>>,
    calc_tt: Option<&'a OrbitalCalcTt<'a>>,
    orbital_ut1_to_tt: Option<&'a dyn Fn(f64) -> Result<f64>>,
    orbital_icrs_at: Option<&'a OrbitalIcrsAt<'a>>,
    orbital_emb_at: Option<&'a OrbitalEmbAt<'a>>,
    orbital_target_lookup: Option<&'a OrbitalTargetLookup<'a>>,
    current_calc_tt: Option<&'a OrbitalCalcTt<'a>>,
    nodaps_target_lookup: Option<&'a OrbitalTargetLookup<'a>>,
    nodaps_ut1_at: Option<&'a dyn Fn(f64) -> Result<f64>>,
    nodaps_calc_ut: Option<&'a NodapsCalcUt<'a>>,
    nodaps_deflect_at: Option<&'a NodapsDeflectorAt<'a>>,
    nodaps_site_state_at: Option<&'a NodapsSiteStateAt<'a>>,
    nodaps_moon_state_at: Option<&'a NodapsMoonStateAt<'a>>,
    nodaps_state_error: Option<&'a dyn Fn(Error) -> Error>,
    smooth_support: std::cell::Cell<Option<NodapsSmoothSupport>>,
}

/// Public TT position feed: scalar Julian Day, body and calculation flags
/// produce the six position components and returned flags, or a typed error.
type OrbitalCalcTt<'a> = dyn Fn(f64, i32, i32) -> Result<(PositionResult, i32)> + 'a;

/// Genuine public lunar-point calculation at scalar UT1 Julian Day, body ID
/// and calculation flags. Returns degrees, AU and per-day spherical rates
/// plus returned flags; public model clocks, guards and errors belong to it.
type NodapsCalcUt<'a> = dyn Fn(f64, i32, i32) -> Result<(PositionResult, i32)> + 'a;

/// Fresh scalar-TT native topocenter-minus-geocenter ICRS state in AU/AU-day.
/// The site is east longitude and latitude in degrees, height in metres.
/// Each call resolves Earth and evaluates Earth+Topos before Earth alone,
/// using independently acquired global UT1. Raw registration/directory and
/// non-finite state errors propagate; no public scalar-TT guard is added.
type NodapsSiteStateAt<'a> = dyn Fn(f64, (f64, f64, f64)) -> Result<([f64; 3], [f64; 3])> + 'a;

/// Fresh scalar-TT Moon-minus-Earth state in true-ecliptic-of-date AU/AU-day.
/// Ordered raw SPK constituents and independent frame clocks belong to the
/// feed; registration, directory, frame and finite-state errors propagate.
type NodapsMoonStateAt<'a> = dyn Fn(f64) -> Result<crate::lunar::EclipticMoonState> + 'a;

/// Raw SSB geometric ICRS acquisition at a scalar TT observable. The provider
/// retains any independently constructed original Time parts internally.
type OrbitalIcrsAt<'a> = dyn Fn(i32, f64) -> Result<crate::positions::GeometricState> + 'a;
/// Direct Earth-Moon-barycenter geometric acquisition at the same TT observable.
type OrbitalEmbAt<'a> = dyn Fn(f64) -> Result<crate::positions::GeometricState> + 'a;
/// Resolve the orbital target alias without evaluating its epoch or center.
type OrbitalTargetLookup<'a> = dyn Fn(i32) -> Result<()> + 'a;
/// Optional raw SSB ICRS deflector position in AU at scalar TT.
type NodapsDeflectorAt<'a> = dyn Fn(i32, f64) -> Result<Option<[f64; 3]>> + 'a;

/// Public orbital-element acquisition error semantics. A missing alias
/// identifies the absent stored support body. A coverage failure retains
/// the failing state's epoch and interval while naming the requested
/// orbital target, including a failure while evaluating its center.
fn orbital_element_state_error(target: i32, error: Error) -> Error {
    match error {
        Error::LebBodyNotFound { body_id } => Error::UnknownBody { body_id },
        Error::LebOutOfRange { jd, start, end, .. } => Error::EphemerisRange {
            requested_jd: jd,
            start_jd: start,
            end_jd: end,
            body: Some(target),
        },
        other => other,
    }
}

impl<'a, P> SealedOrbitalServing<'a, P>
where
    P: StateProvider,
{
    /// Borrow the sealed serving pair and model Delta-T callback.
    ///
    /// `deltat_at(jd)` returns TT-minus-UT1 in **days**, not seconds, at the
    /// supplied scalar Julian day. UT entries evaluate it at UT1 and add the
    /// offset; model terms that require an inverse may evaluate it at TT.
    /// [`Self::with_orbital_ut1_to_tt`] can select a separate clock for orbital
    /// elements without changing node/apse or distance-extrema conversion.
    /// Construction invokes no callback and performs no validation or I/O;
    /// each consuming operation checks its own object, epoch and support domain.
    pub fn new(
        serving: &'a SealedServing<'a, P>,
        deltat_at: &'a dyn Fn(f64) -> Result<f64>,
    ) -> Self {
        Self {
            serving,
            deltat_at,
            ayanamsha: None,
            calc_tt: None,
            orbital_ut1_to_tt: None,
            orbital_icrs_at: None,
            orbital_emb_at: None,
            orbital_target_lookup: None,
            current_calc_tt: None,
            nodaps_target_lookup: None,
            nodaps_ut1_at: None,
            nodaps_calc_ut: None,
            nodaps_deflect_at: None,
            nodaps_site_state_at: None,
            nodaps_moon_state_at: None,
            nodaps_state_error: None,
            smooth_support: std::cell::Cell::new(None),
        }
    }

    /// Attach public TT position acquisition for minor and fictitious bodies.
    ///
    /// The callback receives scalar TT JD, body ID and the exact calculation
    /// flags chosen by the element reducer or current-distance request. It
    /// returns the public six-component position and returned flags. Spherical
    /// reduction requests use degrees, AU and per-day rates; current-distance
    /// requests retain the caller's public output flags. The callback owns
    /// model/provider priority, configuration, frame clocks and typed errors.
    /// Ordinary stored planet/Moon acquisition remains unchanged. Without this
    /// feed, the same calls use the bound stored serving path. No callback
    /// result is cached, and the borrowed feed cannot outlive this serving.
    pub fn with_calc_tt(mut self, calc_tt: &'a OrbitalCalcTt<'a>) -> Self {
        self.calc_tt = Some(calc_tt);
        self
    }

    /// Attach checked raw orbital geometry from an explicitly selected provider.
    ///
    /// `icrs_at(body, tt)` returns a raw SSB-relative ICRS geometric state in AU
    /// and AU/day. `emb_at(tt)` returns the Earth-Moon system barycenter in the
    /// same frame and units. The orbital target resolver runs before evaluation
    /// of the target's center, preserving missing-alias priority. It resolves
    /// Earth to the EMB alias, and every other ordinary target to its source
    /// center-first/barycenter alias. Minor/fictitious reductions instead use
    /// the public calculation feed. Every state must have the requested scalar
    /// TT observable, finite components and an ICRS frame; the shared checked
    /// state access rejects invalid frames and non-finite components.
    ///
    /// The feed owns its source Time representation. In particular, native SPK
    /// public UT reads retain the original UT1 value plus its independently
    /// refined TT fraction, while scalar TT reads preserve their own Time
    /// construction. Raw states enforce their provider directory domains, not
    /// the separate public `calc` scalar-TT envelope. The feed takes no ownership
    /// and cannot outlive this serving. Without it, stored LEB acquisition and
    /// Earth/Moon EMB composition remain unchanged.
    pub fn with_orbital_geometry(
        mut self,
        icrs_at: &'a OrbitalIcrsAt<'a>,
        emb_at: &'a OrbitalEmbAt<'a>,
        target_lookup: &'a OrbitalTargetLookup<'a>,
    ) -> Self {
        self.orbital_icrs_at = Some(icrs_at);
        self.orbital_emb_at = Some(emb_at);
        self.orbital_target_lookup = Some(target_lookup);
        self
    }

    /// Attach the genuine public TT calculation for every extrema current value.
    ///
    /// The callback receives the original scalar TT, requested body and the
    /// source-normalized caller flags with TRUEPOS added. Component 2 returns
    /// after public frame/output shaping, including signed Cartesian z.
    /// This feed owns public scalar-TT guards and configuration/flag errors;
    /// it is distinct from the raw geometry used to fit the osculating ellipse.
    /// Without it, the existing stored-planet and model-only callback policies
    /// remain unchanged. No result or borrowed callback escapes the calculation.
    pub fn with_current_calc_tt(mut self, calc_tt: &'a OrbitalCalcTt<'a>) -> Self {
        self.current_calc_tt = Some(calc_tt);
        self
    }

    /// Attach the native public node/apse acquisition channels.
    ///
    /// Raw body/EMB feeds return finite SSB ICRS states in AU and AU/day at
    /// their scalar TT observable, while retaining the original Time parts.
    /// `target_lookup` resolves geocenter Earth and source planet aliases;
    /// the direct EMB feed owns its separate NAIF 3 alias. `ut1_at` reads
    /// the observation Time's global UT1, retaining original UT1 centrally
    /// and using fresh scalar TT Times for speed sides. `calc_ut` supplies
    /// genuine public spherical lunar point requests; the separately attached
    /// TT feed supplies minor-body states. Model clocks, flags and errors stay
    /// in those public callbacks. `deflector_at` implements the source raw
    /// scalar-Time Sun/Jupiter/Saturn barycenter lookup and optional miss policy.
    ///
    /// These feeds do not apply public calc coverage to raw node/apse geometry.
    /// Native frames/nutation, model selection, point projection and output
    /// shaping remain the shared entry's job. With no feeds, legacy sealed
    /// acquisition and caller-model/matched-epoch clocks stay unchanged.
    /// The private invariant is that these native feeds are attached together;
    /// a lunar public callback therefore always has its Time.ut 1 feed. Every
    /// callback is borrowed for this serving's lifetime; no I/O or ownership
    /// transfer occurs here.
    pub fn with_nodaps_acquisition(
        mut self,
        icrs_at: &'a OrbitalIcrsAt<'a>,
        emb_at: &'a OrbitalEmbAt<'a>,
        target_lookup: &'a OrbitalTargetLookup<'a>,
        ut1_at: &'a dyn Fn(f64) -> Result<f64>,
        calc_ut: &'a NodapsCalcUt<'a>,
        deflector_at: &'a NodapsDeflectorAt<'a>,
    ) -> Self {
        self.orbital_icrs_at = Some(icrs_at);
        self.orbital_emb_at = Some(emb_at);
        self.nodaps_target_lookup = Some(target_lookup);
        self.nodaps_ut1_at = Some(ut1_at);
        self.nodaps_calc_ut = Some(calc_ut);
        self.nodaps_deflect_at = Some(deflector_at);
        self
    }

    /// Attach original observation UT1 and genuine public lunar point reads.
    /// The central UT entry retains its caller UT1 value; speed sides supply
    /// fresh scalar-TT inverse clocks. The two callbacks are installed together
    /// so every lunar `calc_ut` has its corresponding acquired Time observable.
    /// Caller model settings remain attached to the public read. This does not
    /// replace raw body/EMB geometry or change callback-absent serving policy.
    /// Callbacks are borrowed for this synchronous serving lifetime only.
    pub fn with_nodaps_lunar_clock(
        mut self,
        ut1_at: &'a dyn Fn(f64) -> Result<f64>,
        calc_ut: &'a NodapsCalcUt<'a>,
    ) -> Self {
        self.nodaps_ut1_at = Some(ut1_at);
        self.nodaps_calc_ut = Some(calc_ut);
        self
    }

    /// Attach the scalar vector provider's error constructor boundary.
    /// The callback consumes an actual failed raw state and may add owned
    /// constructor facts without changing the native failure. Only node/apse
    /// support and deflector reads use it; ordinary orbital-element policy
    /// remains separate. No error is inferred from output or diagnostic text.
    pub(crate) fn with_nodaps_state_error(mut self, at: &'a dyn Fn(Error) -> Error) -> Self {
        self.nodaps_state_error = Some(at);
        self
    }

    /// Attach native topocentric state acquisition for node/apse output.
    ///
    /// The callback receives fresh scalar TT and a validated geographic site;
    /// it returns topocenter-minus-geocenter position/velocity in ICRS AU and
    /// AU/day, preserving Earth alias/read domains and independent global UT1.
    /// Every nonzero-distance point calls it centrally. Apparent output calls
    /// it again at central and ±0.0001-day epochs, regardless of SPEED; TRUEPOS
    /// and NOABERR skip those extra reads. Zero sentinels skip it; NaN points
    /// still acquire support and retain the output reducer's NaN semantics.
    /// Without this feed the stored site's accepted acquisition is unchanged.
    /// The borrowed callback cannot outlive the serving or transfer ownership.
    pub fn with_nodaps_site_state(mut self, at: &'a NodapsSiteStateAt<'a>) -> Self {
        self.nodaps_site_state_at = Some(at);
        self
    }

    /// Attach the native true-ecliptic Moon state for osculating perigee.
    /// The callback takes fresh scalar TT and returns geometric position and
    /// physical velocity in AU/AU-day. It owns ordered constituent acquisition,
    /// independent TDB/frame clocks and raw error domains; no public calc range
    /// is added. The stored callback-absent lunar acquisition remains unchanged.
    /// This borrowed feed cannot outlive the serving or transfer ownership.
    pub(crate) fn with_nodaps_moon_state(mut self, at: &'a NodapsMoonStateAt<'a>) -> Self {
        self.nodaps_moon_state_at = Some(at);
        self
    }

    /// Select the UT1-to-TT clock for orbital elements only.
    ///
    /// Both input and output are scalar Julian Days. Production public orbital
    /// elements attach the independent global timescale; low-level callers
    /// without this feed retain `UT1 + deltat_at(UT1)`, including matched-epoch
    /// diagnostics. Object validation precedes the conversion, and non-finite
    /// input/output epochs and callback errors refuse the calculation.
    /// [`Self::jd_tt_for_ut`], node/apse UT entries and the distance-extrema UT
    /// extension always retain the bound caller-model Delta-T clock.
    pub fn with_orbital_ut1_to_tt(mut self, ut1_to_tt: &'a dyn Fn(f64) -> Result<f64>) -> Self {
        self.orbital_ut1_to_tt = Some(ut1_to_tt);
        self
    }

    /// Attach the sealed ayanamsha serving: enables the plain-`FLG_SIDEREAL` mean-ayanamsha arm. The
    /// serving owns selection and configured USER_UT/Valens defining epochs.
    /// Native node/apse observation ayanamsha instead uses global TT acquired
    /// from Time.ut 1; callback-absent stored calls retain this serving's
    /// configured observation conversion. No I/O or globals.
    pub fn with_ayanamsha(mut self, ayanamsha: SealedAyanamshaServing<'a>) -> Self {
        self.ayanamsha = Some(ayanamsha);
        self
    }

    /// Bound sealed serving pair.
    pub fn serving(&self) -> &SealedServing<'a, P> {
        self.serving
    }

    /// Convert a finite UT1 Julian day to TT using `UT1 + deltat_at(UT1)`.
    ///
    /// The callback returns days; callback errors and nonfinite input/offset
    /// values are errors. This uses the bound model callback even when an
    /// independent orbital-element clock has been attached.
    pub fn jd_tt_for_ut(&self, jd_ut: f64) -> Result<f64> {
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let dt = (self.deltat_at)(jd_ut)?;
        if !dt.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt });
        }
        Ok(jd_ut + dt)
    }

    /// Intersect an actually consumed polynomial's domain in the outer TT clock.
    fn note_smooth_interval(&self, query_tt: f64, interval: Option<(f64, f64)>) {
        if let (Some(mut support), Some((start, end))) = (self.smooth_support.get(), interval) {
            let offset = query_tt - support.epoch;
            support.start = support.start.max(start - offset);
            support.end = support.end.min(end - offset);
            self.smooth_support.set(Some(support));
        }
    }

    fn note_body_support(&self, body: i32, query_tt: f64) -> Result<()> {
        if self.smooth_support.get().is_some() {
            self.note_smooth_interval(
                query_tt,
                self.serving.provider().smooth_interval_tt(body, query_tt)?,
            );
        }
        Ok(())
    }

    fn note_nutation_support(&self, query_tt: f64) -> Result<()> {
        if self.smooth_support.get().is_some() {
            self.note_smooth_interval(
                query_tt,
                self.serving
                    .provider()
                    .nutation_smooth_interval_tt(query_tt)?,
            );
        }
        Ok(())
    }

    /// Geometric ICRS position in AU and velocity in AU/day. An attached raw
    /// orbital callback owns acquisition and its provider directory domain;
    /// the legacy stored-state branch first enforces the serving channel gate.
    /// Both branches require the requested body, identical scalar TT epoch,
    /// ICRS frame and finite components. Acquisition and invariant errors
    /// propagate before orbital reduction.
    fn icrs_state(&self, body_id: i32, jd_tt: f64) -> Result<([f64; 3], [f64; 3])> {
        self.note_body_support(body_id, jd_tt)?;
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let state = if let Some(icrs_at) = self.orbital_icrs_at {
            icrs_at(body_id, jd_tt)?
        } else {
            serving_channel_gate(self.serving.provider(), body_id)?;
            self.serving.provider().geometric_state(body_id, jd_tt)?
        };
        self.checked_orbital_icrs_state(body_id, jd_tt, state)
    }

    /// Same-epoch geometric displacement and velocity in ICRS AU/AU-day.
    ///
    /// Stored providers can subtract aligned polynomial coefficients before
    /// evaluation, avoiding cancellation between rounded barycentric states.
    /// Native acquisition callbacks retain their explicit clock and alias
    /// policy. Both paths validate the requested states before conic reduction.
    fn relative_icrs_state(
        &self,
        target: i32,
        center: i32,
        jd_tt: f64,
    ) -> Result<([f64; 3], [f64; 3])> {
        self.note_body_support(target, jd_tt)?;
        self.note_body_support(center, jd_tt)?;
        if self.orbital_icrs_at.is_some() {
            let (target_p, target_v) = self.icrs_state(target, jd_tt)?;
            let (center_p, center_v) = self.icrs_state(center, jd_tt)?;
            return Ok((
                std::array::from_fn(|axis| target_p[axis] - center_p[axis]),
                std::array::from_fn(|axis| target_v[axis] - center_v[axis]),
            ));
        }
        serving_channel_gate(self.serving.provider(), target)?;
        serving_channel_gate(self.serving.provider(), center)?;
        let state = self
            .serving
            .provider()
            .relative_geometric_state(target, center, jd_tt)?;
        if state.body_id() != target || state.center_body_id() != center {
            return Err(Error::InvalidOrbitalObject {
                entry: "orbital_geometry",
                detail: "relative state identifies different bodies".to_string(),
            });
        }
        if state.jd_tt().to_bits() != jd_tt.to_bits() {
            return Err(Error::MathDomain {
                value: state.jd_tt(),
            });
        }
        let (position, velocity) = (state.pos(), state.vel());
        for value in position.into_iter().chain(velocity) {
            if !value.is_finite() {
                return Err(Error::NonFiniteArgument { value });
            }
        }
        Ok((position, velocity))
    }

    /// Validate frame, scalar epoch identity and finite raw components before
    /// any orbital arithmetic. Typed provider acquisition errors propagate
    /// before this invariant check; no apparent-place corrections occur here.
    fn checked_orbital_icrs_state(
        &self,
        body_id: i32,
        jd_tt: f64,
        state: crate::positions::GeometricState,
    ) -> Result<([f64; 3], [f64; 3])> {
        if state.body_id() != body_id {
            return Err(Error::InvalidOrbitalObject {
                entry: "orbital_geometry",
                detail: format!(
                    "state for body {} does not identify requested body {body_id}",
                    state.body_id()
                ),
            });
        }
        if state.jd_tt().to_bits() != jd_tt.to_bits() {
            return Err(Error::MathDomain {
                value: state.jd_tt(),
            });
        }
        if !state.frame().is_icrs() {
            return Err(Error::SealedChannelNotServed {
                body_id,
                coord_type: self
                    .serving
                    .provider()
                    .serving_coord_type(body_id)
                    .unwrap_or(0),
            });
        }
        let (p, v) = (state.pos(), state.vel());
        for c in p.into_iter().chain(v) {
            if !c.is_finite() {
                return Err(Error::NonFiniteArgument { value: c });
            }
        }
        Ok((p, v))
    }

    /// Osculating elements at a TT epoch. `flags` carries `FLG_BARYCTR` /
    /// `FLG_ORBEL_AA` (bit 32768 reads as `FLG_ORBEL_AA` here, since
    /// topocentric places have no meaning for elements). Fictitious bodies
    /// ignore BARYCTR. White Moon/Waldemath additionally ignore ORBEL_AA and
    /// retain Earth-centered J2000 geometry; Proserpina retains the circular
    /// mean-of-date geometry. Models are obtained from the attached public
    /// calculation feed, or retain the stored-channel refusal without it.
    /// Ordinary planet/Moon reductions refer directions to the J2000 ecliptic
    /// and equinox; their center and effective mass follow the selected body
    /// and flags. The model-body exceptions above retain their stated frames.
    /// The 50 returned slots and their exact units are documented by
    /// [`orbital_elements_from_ecliptic_state`]; slots 17–49 are positive zeros.
    /// Invalid orbital objects, nonfinite epochs and consumed state/frame or
    /// coverage failures return typed errors, without a substitute element set.
    pub fn orbital_elements_tt(&self, body_id: i32, jd_tt: f64, flags: i32) -> Result<[f64; 50]> {
        validate_orbital_object("get_orbital_elements", body_id, false)?;
        self.elements_inner(body_id, jd_tt, flags)
    }

    /// Osculating elements at a UT1 JD. An attached orbital timescale converts
    /// to TT; otherwise the bound caller-model Delta-T supplies TT. The public
    /// assembly attaches the global timescale, while matched-epoch low-level
    /// diagnostics retain their explicit caller-model conversion.
    /// The 50 slots and units follow [`orbital_elements_from_ecliptic_state`].
    pub fn orbital_elements_ut(&self, body_id: i32, jd_ut: f64, flags: i32) -> Result<[f64; 50]> {
        validate_orbital_object("get_orbital_elements_ut", body_id, false)?;
        let jd_tt = if let Some(to_tt) = self.orbital_ut1_to_tt {
            if !jd_ut.is_finite() {
                return Err(Error::NonFiniteArgument { value: jd_ut });
            }
            let tt = to_tt(jd_ut)?;
            if !tt.is_finite() {
                return Err(Error::NonFiniteArgument { value: tt });
            }
            tt
        } else {
            self.jd_tt_for_ut(jd_ut)?
        };
        self.elements_inner(body_id, jd_tt, flags)
    }

    /// Shared element reduction behind both entries.
    fn elements_inner(&self, body_id: i32, jd_tt: f64, flags: i32) -> Result<[f64; 50]> {
        // Both public entries reject the Sun before reaching this reduction:
        // the Sun has no orbital-element tuple.
        if is_fict_geocentric_body(body_id) {
            return self.calc_state_geocentric_elements(body_id, jd_tt, flags);
        }
        if is_unserved_model_body(body_id) {
            // Minor/fictitious reductions read the geometric public position.
            // A bound callback supplies runtime models; without one the stored
            // serving's typed channel refusal still propagates.
            return self.calc_state_elements(body_id, jd_tt, flags);
        }
        if !is_nodaps_planet(body_id) {
            return Err(Error::InvalidOrbitalObject {
                entry: "get_orbital_elements",
                detail: format!("error in get_orbital_elements(): object {body_id} not valid"),
            });
        }
        // Target aliases are resolved before their center is evaluated.
        // The Earth orbital alias requires both Earth and Moon channels.
        if let Some(target_lookup) = self.orbital_target_lookup {
            target_lookup(body_id)?;
        } else {
            let required_targets: &[i32] = if body_id == constants::EARTH {
                &[constants::EARTH, constants::MOON]
            } else {
                std::slice::from_ref(&body_id)
            };
            for &target in required_targets {
                crate::positions::serving_channel_gate(self.serving.provider(), target)
                    .map_err(|error| orbital_element_state_error(body_id, error))?;
            }
        }
        if body_id == constants::MOON {
            return self.moon_elements_tt(jd_tt);
        }
        let (r_helio, v_helio) = self.helio_icrs_state(body_id, jd_tt)?;
        let gm = if (flags & constants::FLG_ORBEL_AA) != 0 {
            match orbel_aa_cumulative_through_self(body_id) {
                Some(excess) => GM_SUN * (1.0 + excess),
                None => planet_system_gm(body_id).unwrap_or(GM_SUN),
            }
        } else {
            planet_system_gm(body_id).unwrap_or(GM_SUN)
        };
        let m = icrs_to_ecl_j2000_matrix();
        let pos = mat3_vec3(m, r_helio)?;
        let vel = mat3_vec3(m, v_helio)?;
        let (mut x, mut y, mut z, mut vx, mut vy, mut vz) =
            (pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]);
        // Barycentric re-reference for trans-jovian orbits only: the
        // classifier runs on the default solar-GM semi-major axis, and
        // the SSB state is already SSB-relative, so only the J2000
        // frame rotation applies.
        if (flags & constants::FLG_BARYCTR) != 0
            && osculating_semi_major_axis(x, y, z, vx, vy, vz, gm) > TRANS_JOVIAN_A_AU
        {
            let (tp, tv) = self.icrs_state(body_id, jd_tt)?;
            let shifted = mat3_vec3(m, tp)?;
            let shifted_v = mat3_vec3(m, tv)?;
            x = shifted[0];
            y = shifted[1];
            z = shifted[2];
            vx = shifted_v[0];
            vy = shifted_v[1];
            vz = shifted_v[2];
        }
        Ok(orbital_elements_from_ecliptic_state(
            [x, y, z, vx, vy, vz],
            gm,
            body_id,
            jd_tt,
            Some(GM_SUN),
        ))
    }

    /// Heliocentric ICRS state in AU and AU/day (target minus Sun).
    /// The Earth orbital target is the Earth-Moon barycenter; all other
    /// targets use their own stored ICRS channel. Support-channel errors
    /// propagate before the element reduction.
    fn helio_icrs_state(&self, target: i32, jd_tt: f64) -> Result<([f64; 3], [f64; 3])> {
        let (sp, sv) = self
            .icrs_state(constants::SUN, jd_tt)
            .map_err(|error| orbital_element_state_error(target, error))?;
        let (tp, tv) = if target == constants::EARTH {
            self.emb_icrs_state(jd_tt)
        } else {
            self.icrs_state(target, jd_tt)
        }
        .map_err(|error| orbital_element_state_error(target, error))?;
        Ok((
            [tp[0] - sp[0], tp[1] - sp[1], tp[2] - sp[2]],
            [tv[0] - sv[0], tv[1] - sv[1], tv[2] - sv[2]],
        ))
    }

    /// Geocentric lunar elements at a TT epoch.
    fn moon_elements_tt(&self, jd_tt: f64) -> Result<[f64; 50]> {
        // Direct native callers (including extrema) preserve the source's
        // target alias before its Earth center. The legacy stored path keeps
        // its existing channel policy; elements_inner already gates its target.
        if let Some(target_lookup) = self.orbital_target_lookup {
            target_lookup(constants::MOON)?;
        }
        let (r, v) = self
            .relative_icrs_state(constants::MOON, constants::EARTH, jd_tt)
            .map_err(|error| orbital_element_state_error(constants::MOON, error))?;
        let m = icrs_to_ecl_j2000_matrix();
        let pos = mat3_vec3(m, r)?;
        let vel = mat3_vec3(m, v)?;
        Ok(orbital_elements_from_ecliptic_state(
            [pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]],
            moon_system_gm(),
            constants::MOON,
            jd_tt,
            Some(moon_central_gm()),
        ))
    }

    /// Acquire a minor/fictitious body's public scalar-TT position, or retain
    /// the stored serving policy when no public calculation feed is attached.
    /// Callback flags and error categories are passed through without rewriting.
    fn model_calc_tt(&self, jd_tt: f64, body_id: i32, flags: i32) -> Result<(PositionResult, i32)> {
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        if let Some(calc_tt) = self.calc_tt {
            calc_tt(jd_tt, body_id, flags)
        } else {
            serving_compute_tt(self.serving, jd_tt, body_id, flags)
        }
    }

    /// Reduce a geometric heliocentric public position with solar GM.
    /// Minor bodies use the J2000 ecliptic and may re-reference trans-jovian
    /// states to the SSB. Fictitious bodies ignore BARYCTR; Proserpina uses its
    /// native mean ecliptic of date to retain a self-consistent circular speed.
    /// ORBEL_AA augments the geometry GM while periods retain pure solar GM.
    /// Positions are degrees/AU with per-day rates; degenerate distances and
    /// public acquisition errors refuse before reduction.
    fn calc_state_elements(&self, body_id: i32, jd_tt: f64, flags: i32) -> Result<[f64; 50]> {
        let fictitious = (constants::CUPIDO..=constants::WALDEMATH).contains(&body_id);
        let helio_flags = constants::FLG_HELCTR
            | if body_id == constants::PROSERPINA {
                constants::FLG_NONUT
            } else {
                constants::FLG_J2000
            }
            | constants::FLG_TRUEPOS
            | constants::FLG_SPEED;
        let (pos, _) = self.model_calc_tt(jd_tt, body_id, helio_flags)?;
        if !pos[2].is_finite() || pos[2] <= 0.0 {
            return Err(Error::NonFiniteArgument { value: pos[2] });
        }
        let cart = ecliptic_spherical_to_cartesian_state(pos);
        let (mut x, mut y, mut z, mut vx, mut vy, mut vz) =
            (cart[0], cart[1], cart[2], cart[3], cart[4], cart[5]);
        let gm = GM_SUN;
        let mut gm_slots: Option<f64> = None;
        if !fictitious
            && (flags & constants::FLG_BARYCTR) != 0
            && osculating_semi_major_axis(x, y, z, vx, vy, vz, gm) > TRANS_JOVIAN_A_AU
        {
            // SSB re-reference by adding the Sun's barycentric
            // J2000-ecliptic state; the sealed Sun state
            // is already SSB-relative, so only the frame rotation
            // applies.
            let (sp, sv) = self.icrs_state(constants::SUN, jd_tt)?;
            let m = icrs_to_ecl_j2000_matrix();
            let s = mat3_vec3(m, sp)?;
            let vs = mat3_vec3(m, sv)?;
            x += s[0];
            y += s[1];
            z += s[2];
            vx += vs[0];
            vy += vs[1];
            vz += vs[2];
        }
        let mut gm_geo = gm;
        if (flags & constants::FLG_ORBEL_AA) != 0 {
            let a_default = osculating_semi_major_axis(x, y, z, vx, vy, vz, gm);
            gm_slots = Some(gm);
            gm_geo = GM_SUN * (1.0 + orbel_aa_interior_mass_sum(a_default));
        }
        Ok(orbital_elements_from_ecliptic_state(
            [x, y, z, vx, vy, vz],
            gm_geo,
            body_id,
            jd_tt,
            gm_slots,
        ))
    }

    /// Ascending node, descending node, periapsis and apoapsis at a TT Julian day.
    ///
    /// Each point is `[longitude, latitude, distance, longitude_rate,
    /// latitude_rate, radial_rate]`, in degrees, AU and corresponding units
    /// per TT day. Without `FLG_SPEED`, rates are zero. `FLG_EQUATORIAL` selects
    /// RA/declination in degrees; XYZ selects AU/AU-day and takes precedence
    /// over RADIANS, which otherwise changes angular slots to radians.
    ///
    /// [`prefer_mean`] defines the method's low-byte mean/osculating choice;
    /// bodies without a mean-element table use osculating geometry.
    /// `NODBIT_OSCU_BAR` enables the trans-Jovian barycentric fit, while
    /// `NODBIT_FOPOINT` replaces the apoapsis with the second focal point.
    /// These are geometric constructions, not an integrated orbit prediction.
    ///
    /// The ordinary working frame is true ecliptic/equinox of date using native
    /// IAU 2006/2000A nutation; NONUT selects the mean frame. J2000 and sidereal
    /// flags apply their declared frame reductions. SIDEREAL requires an attached
    /// ayanamsha serving, including fixed-epoch and active-plane projections.
    /// Centered views suppress topocentric parallax; otherwise TOPOCTR requires
    /// [`Self::nod_aps_topo_tt`] and its explicit site.
    ///
    /// SPEED normally uses fourth-order differences within the common consumed
    /// polynomial support: centered where possible, one-sided near an edge, with
    /// a reduced step when necessary. The mean-Moon arm can supply analytic rates.
    /// Insufficient representable support is an error, not an extrapolated rate.
    /// Undefined Sun/Earth nodes use zero sentinels; interpolated lunar apsides
    /// use [`intp_undefined_points`] NaNs rather than physical node/apse positions.
    /// Unsupported bodies/flags, nonfinite inputs and consumed provider, clock,
    /// frame or coverage failures return typed errors.
    pub fn nod_aps_tt(
        &self,
        body_id: i32,
        jd_tt: f64,
        method: i32,
        flags: i32,
    ) -> Result<[PositionResult; 4]> {
        // Parallax has no effect on centered heliocentric/barycentric
        // points, so those requests do not require an observer.
        let flags = if (flags & (constants::FLG_HELCTR | constants::FLG_BARYCTR)) != 0 {
            flags & !constants::FLG_TOPOCTR
        } else {
            flags
        };
        if self.nodaps_projection_route(flags).is_some() {
            return self.nod_aps_projected_tt(body_id, jd_tt, method, flags);
        }
        if self.nodaps_target_lookup.is_some() && (flags & constants::FLG_TOPOCTR) != 0 {
            // Source formatting checks the observer only after computing the
            // unshaped points; zero sentinels do not require a site. Projected
            // sub-requests re-enter here before their rigid rotation.
            let pts =
                self.nod_aps_inner_tt(body_id, jd_tt, method, flags & !constants::FLG_TOPOCTR)?;
            if pts.iter().any(|point| point[2] != 0.0) {
                return Err(Error::ConfigurationMissing {
                    missing_config: "observer_location",
                    suggestion: "call set_topo(lon, lat, alt) first",
                });
            }
            let pts = self.apply_nodaps_equatorial(pts, jd_tt, flags)?;
            return Ok(shape_nodaps_output(pts, flags));
        }
        self.nod_aps_plain_tt(body_id, jd_tt, method, flags)
    }

    /// Unprojected node/apse 4-tuple at a TT epoch: the inner spherical
    /// computation, the entry-level equatorial rotation, then the
    /// output-flag shaping. The
    /// projected arm re-enters through [`nod_aps_tt`](Self::nod_aps_tt)
    /// with `FLG_SIDEREAL` stripped, so it always lands here.
    fn nod_aps_plain_tt(
        &self,
        body_id: i32,
        jd_tt: f64,
        method: i32,
        flags: i32,
    ) -> Result<[PositionResult; 4]> {
        let pts = self.nod_aps_inner_tt(body_id, jd_tt, method, flags)?;
        let pts = self.apply_nodaps_equatorial(pts, jd_tt, flags)?;
        Ok(shape_nodaps_output(pts, flags))
    }

    /// Ecliptic-to-equatorial rotation at the entry. The
    /// rotation covers the speeds with the state because the stencil
    /// settled in the ecliptic frame, like the calculation order. Obliquity
    /// selection drops the sidereal bit first: J2000 legs use the IAU 2006 J2000 constant,
    /// `FLG_NONUT` the mean obliquity, otherwise the true obliquity
    /// (mean plus the native nutation-in-obliquity value, the same
    /// `erfa.nut06a` spelling the working frame evaluates).
    /// Non-finite points pass through.
    fn apply_nodaps_equatorial(
        &self,
        pts: [PositionResult; 4],
        jd_tt: f64,
        flags: i32,
    ) -> Result<[PositionResult; 4]> {
        if (flags & constants::FLG_EQUATORIAL) == 0 {
            return Ok(pts);
        }
        let eq = flags & !constants::FLG_SIDEREAL;
        let eps_deg = if (eq & constants::FLG_J2000) != 0 {
            // IAU 2006 mean obliquity at J2000.
            84381.406 / 3600.0
        } else {
            let mean = vondrak_mean_obliquity_deg(jd_tt)?;
            if (eq & constants::FLG_NONUT) != 0 {
                mean
            } else {
                mean + crate::erfa::nut06a(jd_tt)?.1.to_degrees()
            }
        };
        let mut out = pts;
        for pt in out.iter_mut() {
            if !(pt[0].is_finite() && pt[1].is_finite()) {
                continue;
            }
            let ([ra, dec, dist], [dra, ddec, ddist]) =
                crate::utils::cotrans_sp(pt[0], pt[1], pt[2], pt[3], pt[4], pt[5], -eps_deg)?;
            *pt = [ra, dec, dist, dra, ddec, ddist];
        }
        Ok(out)
    }

    /// ICRS observer offset mapped into the node/apse output frame at one
    /// epoch: the
    /// Vondrák precession-nutation matrix plus the matching obliquity
    /// (mean frame under `FLG_NONUT`, native IAU 2006/2000 A pair
    /// otherwise, evaluated per epoch like the working frame), then the
    /// J2000 vector rotation and the true-ayanamsha rotation for the
    /// sidereal frame. Linear, so position and velocity vectors
    /// transform identically.
    fn nodaps_output_frame_vector(
        &self,
        v: [f64; 3],
        epoch_tt: f64,
        flags: i32,
    ) -> Result<[f64; 3]> {
        for c in v {
            if !c.is_finite() {
                return Err(Error::NonFiniteArgument { value: c });
            }
        }
        let nut = if (flags & constants::FLG_NONUT) != 0 {
            None
        } else {
            Some(crate::erfa::nut06a(epoch_tt)?)
        };
        let frame = NodapsFrame::for_date(epoch_tt, nut)?;
        let mut e = frame.icrs_to_ecliptic(v);
        if (flags & constants::FLG_J2000) != 0 {
            e = nodaps_vector_to_j2000(e, epoch_tt)?;
        }
        if (flags & constants::FLG_SIDEREAL) != 0 {
            // passes `flags & ~FLG_J2000`
            // (true ayanamsha when neither `FLG_NONUT` nor `FLG_J2000`
            // is set), unlike the body arm below, which forces the
            // mean value. The nutation-in-longitude term is ~arcsec
            // scale on the ~4 e-5 AU observer, i.e. ~1 e-8 deg output —
            // discriminating re-proof of the first-try mean-arm miss.
            let aya_deg = self.nodaps_true_ayanamsha_deg(epoch_tt, flags)?;
            let (sa, ca) = aya_deg.to_radians().sin_cos();
            e = [e[0] * ca + e[1] * sa, -e[0] * sa + e[1] * ca, e[2]];
        }
        Ok(e)
    }

    /// ICRS observer state (position offset and rotation velocity,
    /// AU and AU/day) at one TT epoch under its UT1 partner: the native
    /// `c2t06a` terrestrial-rotation series with zero polar motion (the
    /// served production assembly) through the IERS 2010 site ellipsoid.
    /// This observer policy is distinct from the WGS84 fast-calc channel;
    /// the two ellipsoids must not be interchanged in this reduction.
    fn nodaps_observer_icrs(
        tt: f64,
        ut1: f64,
        site: (f64, f64, f64),
    ) -> Result<([f64; 3], [f64; 3])> {
        if !tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: tt });
        }
        if !ut1.is_finite() {
            return Err(Error::NonFiniteArgument { value: ut1 });
        }
        let c2t = crate::erfa::c2t06a(tt, ut1, 0.0, 0.0)?;
        Ok(topocentric_offset(site.0, site.1, site.2, c2t, SiteEllipsoid::Iers2010)?.pair())
    }

    /// Topocentric parallax for one finished geocentric set: the observer state in the
    /// output frame (physical velocity plus the frame-derivative term
    /// over the ±0.001 d stencil on the fixed central ICRS offset),
    /// subtracted in the ecliptic, with the diurnal-aberration stencil
    /// (±0.0001 d, skipped under `FLG_TRUEPOS`/`FLG_NOABERR` like the
    /// the model) carried into all three speed channels, then the
    /// spherical round-trip with the model guards.
    ///
    /// Centered requests pass through; zero-distance sentinels
    /// (the Sun's undefined node slots) pass through per point, like
    /// the model; non-finite position slots pass through like the
    /// entry equatorial arm (reachable NaN points already carry the
    /// matching zero/NaN speeds, so the values agree with the model
    /// NaN propagation). A zero topocentric distance returns the input
    /// geocentric point, like the model. Speeds are physical rates
    /// only under `FLG_SPEED`, else zeroed. These precomputed state/stencil
    /// semantics describe the callback-absent stored channel. An attached
    /// native site channel instead performs per-point raw support acquisition
    /// and preserves conditional reads and NaN-point formatting through
    /// [`Self::apply_native_nodaps_topocentric`].
    fn apply_nodaps_topocentric(
        &self,
        pts: [PositionResult; 4],
        jd_tt: f64,
        jd_ut1: f64,
        flags: i32,
        site: (f64, f64, f64),
    ) -> Result<[PositionResult; 4]> {
        if (flags & (constants::FLG_HELCTR | constants::FLG_BARYCTR)) != 0 {
            return Ok(pts);
        }
        if let Some(site_state_at) = self.nodaps_site_state_at {
            return self.apply_native_nodaps_topocentric(pts, jd_tt, flags, site, site_state_at);
        }
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        if !jd_ut1.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut1 });
        }
        for v in [site.0, site.1, site.2] {
            if !v.is_finite() {
                return Err(Error::NonFiniteArgument { value: v });
            }
        }
        // Terrestrial-rotation clock: the UT1 partner travels with the
        // stencil (frozen TT-UT1 across ±0.001 d, sub-nanosecond).
        let ut1_at = |tt: f64| jd_ut1 + (tt - jd_tt);
        let (op_icrs, ov_icrs) = Self::nodaps_observer_icrs(jd_tt, jd_ut1, site)?;
        let obs_pos = self.nodaps_output_frame_vector(op_icrs, jd_tt, flags)?;
        let obs_vel_physical = self.nodaps_output_frame_vector(ov_icrs, jd_tt, flags)?;
        // Frame-derivative term on the fixed central ICRS offset.
        let frame_step = 0.001;
        let frame_before = self.nodaps_output_frame_vector(op_icrs, jd_tt - frame_step, flags)?;
        let frame_after = self.nodaps_output_frame_vector(op_icrs, jd_tt + frame_step, flags)?;
        let obs_vel = [
            obs_vel_physical[0] + (frame_after[0] - frame_before[0]) / (2.0 * frame_step),
            obs_vel_physical[1] + (frame_after[1] - frame_before[1]) / (2.0 * frame_step),
            obs_vel_physical[2] + (frame_after[2] - frame_before[2]) / (2.0 * frame_step),
        ];
        // Framed relative (diurnal) velocity at the aberration stencil
        // epochs for the speed-carrying correction below.
        let aber_step = 0.0001;
        let (_, rv_before_icrs) =
            Self::nodaps_observer_icrs(jd_tt - aber_step, ut1_at(jd_tt - aber_step), site)?;
        let (_, rv_after_icrs) =
            Self::nodaps_observer_icrs(jd_tt + aber_step, ut1_at(jd_tt + aber_step), site)?;
        let rv_now_icrs = ov_icrs;
        let rel_now = self.nodaps_output_frame_vector(rv_now_icrs, jd_tt, flags)?;
        let rel_before =
            self.nodaps_output_frame_vector(rv_before_icrs, jd_tt - aber_step, flags)?;
        let rel_after = self.nodaps_output_frame_vector(rv_after_icrs, jd_tt + aber_step, flags)?;
        let want_aberr = (flags & (constants::FLG_TRUEPOS | constants::FLG_NOABERR)) == 0;
        let want_speed = (flags & constants::FLG_SPEED) != 0;
        let mut out = pts;
        for pt in out.iter_mut() {
            let [lon, lat, dist, dlon, dlat, ddist] = *pt;
            if !(lon.is_finite() && lat.is_finite() && dist.is_finite()) {
                continue;
            }
            if dist == 0.0 {
                continue;
            }
            let lon_rad = lon.to_radians();
            let lat_rad = lat.to_radians();
            let dlon_rad = dlon.to_radians();
            let dlat_rad = dlat.to_radians();
            let cl = lon_rad.cos();
            let sl = lon_rad.sin();
            let cb = lat_rad.cos();
            let sb = lat_rad.sin();
            let body_pos = [dist * cb * cl, dist * cb * sl, dist * sb];
            let body_vel = [
                ddist * cb * cl - dist * sb * dlat_rad * cl - dist * cb * sl * dlon_rad,
                ddist * cb * sl - dist * sb * dlat_rad * sl + dist * cb * cl * dlon_rad,
                ddist * sb + dist * cb * dlat_rad,
            ];
            let topo_pos = [
                body_pos[0] - obs_pos[0],
                body_pos[1] - obs_pos[1],
                body_pos[2] - obs_pos[2],
            ];
            let topo_vel = [
                body_vel[0] - obs_vel[0],
                body_vel[1] - obs_vel[1],
                body_vel[2] - obs_vel[2],
            ];
            let (tx, ty, tz, tvx, tvy, tvz) = if want_aberr {
                let aberrated = |pos: [f64; 3], rel: [f64; 3]| apply_aberration(pos, rel, 0.0);
                let p0 = aberrated(topo_pos, rel_now)?;
                let p_before = aberrated(
                    [
                        topo_pos[0] - aber_step * topo_vel[0],
                        topo_pos[1] - aber_step * topo_vel[1],
                        topo_pos[2] - aber_step * topo_vel[2],
                    ],
                    rel_before,
                )?;
                let p_after = aberrated(
                    [
                        topo_pos[0] + aber_step * topo_vel[0],
                        topo_pos[1] + aber_step * topo_vel[1],
                        topo_pos[2] + aber_step * topo_vel[2],
                    ],
                    rel_after,
                )?;
                (
                    p0[0],
                    p0[1],
                    p0[2],
                    (p_after[0] - p_before[0]) / (2.0 * aber_step),
                    (p_after[1] - p_before[1]) / (2.0 * aber_step),
                    (p_after[2] - p_before[2]) / (2.0 * aber_step),
                )
            } else {
                (
                    topo_pos[0],
                    topo_pos[1],
                    topo_pos[2],
                    topo_vel[0],
                    topo_vel[1],
                    topo_vel[2],
                )
            };
            let distance = (tx * tx + ty * ty + tz * tz).sqrt();
            if distance == 0.0 {
                continue;
            }
            let distance_xy_sq = tx * tx + ty * ty;
            let distance_xy = distance_xy_sq.sqrt();
            let radial_speed = (tx * tvx + ty * tvy + tz * tvz) / distance;
            let longitude_speed = if distance_xy_sq > 0.0 {
                ((tx * tvy - ty * tvx) / distance_xy_sq).to_degrees()
            } else {
                0.0
            };
            let latitude_speed = if distance_xy > 0.0 {
                ((tvz - tz * radial_speed / distance) / distance_xy).to_degrees()
            } else {
                0.0
            };
            let (dlon_out, dlat_out, ddist_out) = if want_speed {
                (longitude_speed, latitude_speed, radial_speed)
            } else {
                (0.0, 0.0, 0.0)
            };
            *pt = [
                ty.atan2(tx).to_degrees().rem_euclid(360.0),
                (tz / distance).clamp(-1.0, 1.0).asin().to_degrees(),
                distance,
                dlon_out,
                dlat_out,
                ddist_out,
            ];
        }
        Ok(out)
    }

    /// Native per-point topocentric output with fresh scalar-Time support.
    ///
    /// Only exact zero distance bypasses acquisition. Every other point,
    /// including undefined NaN apsides, acquires the central observer state
    /// before transforming position/physical velocity and the fixed-vector
    /// ±0.001-day frame derivative. Apparent output then reacquires central,
    /// negative and positive ±0.0001-day site states in that order; it uses
    /// the already aberrated central position for the side evaluations.
    /// TRUEPOS/NOABERR omit this branch independently of SPEED. The callback
    /// owns raw Earth errors and each independent global UT1 clock. Result
    /// values are degrees, AU and per-day rates, with source zero/NaN guards.
    fn apply_native_nodaps_topocentric(
        &self,
        pts: [PositionResult; 4],
        jd_tt: f64,
        flags: i32,
        site: (f64, f64, f64),
        site_state_at: &NodapsSiteStateAt<'_>,
    ) -> Result<[PositionResult; 4]> {
        let want_aberr = (flags & (constants::FLG_TRUEPOS | constants::FLG_NOABERR)) == 0;
        let want_speed = (flags & constants::FLG_SPEED) != 0;
        let mut out = pts;
        for point in &mut out {
            let [lon, lat, dist, dlon, dlat, ddist] = *point;
            if dist == 0.0 {
                continue;
            }
            let (op_icrs, ov_icrs) = site_state_at(jd_tt, site)?;
            let obs_pos = self.nodaps_output_frame_vector(op_icrs, jd_tt, flags)?;
            let physical_vel = self.nodaps_output_frame_vector(ov_icrs, jd_tt, flags)?;
            let frame_step = 0.001;
            let before = self.nodaps_output_frame_vector(op_icrs, jd_tt - frame_step, flags)?;
            let after = self.nodaps_output_frame_vector(op_icrs, jd_tt + frame_step, flags)?;
            let obs_vel = std::array::from_fn::<_, 3, _>(|i| {
                physical_vel[i] + (after[i] - before[i]) / (2.0 * frame_step)
            });
            let lon_rad = lon.to_radians();
            let lat_rad = lat.to_radians();
            let dlon_rad = dlon.to_radians();
            let dlat_rad = dlat.to_radians();
            let cl = lon_rad.cos();
            let sl = lon_rad.sin();
            let cb = lat_rad.cos();
            let sb = lat_rad.sin();
            let body_pos = [dist * cb * cl, dist * cb * sl, dist * sb];
            let body_vel = [
                ddist * cb * cl - dist * sb * dlat_rad * cl - dist * cb * sl * dlon_rad,
                ddist * cb * sl - dist * sb * dlat_rad * sl + dist * cb * cl * dlon_rad,
                ddist * sb + dist * cb * dlat_rad,
            ];
            let mut topo_pos = std::array::from_fn::<_, 3, _>(|i| body_pos[i] - obs_pos[i]);
            let mut topo_vel = std::array::from_fn::<_, 3, _>(|i| body_vel[i] - obs_vel[i]);
            if want_aberr {
                let aberrated_at = |base: [f64; 3], offset: f64| -> Result<[f64; 3]> {
                    let position = std::array::from_fn(|i| base[i] + offset * topo_vel[i]);
                    let side_tt = jd_tt + offset;
                    let (_, relative_icrs) = site_state_at(side_tt, site)?;
                    let relative =
                        self.nodaps_output_frame_vector(relative_icrs, side_tt, flags)?;
                    Ok(nodaps_diurnal_aberration(position, relative))
                };
                topo_pos = aberrated_at(topo_pos, 0.0)?;
                let step = 0.0001;
                let before = aberrated_at(topo_pos, -step)?;
                let after = aberrated_at(topo_pos, step)?;
                topo_vel = std::array::from_fn(|i| (after[i] - before[i]) / (2.0 * step));
            }
            let [tx, ty, tz] = topo_pos;
            let [tvx, tvy, tvz] = topo_vel;
            let distance_xy_sq = tx * tx + ty * ty;
            let distance = (distance_xy_sq + tz * tz).sqrt();
            if distance == 0.0 {
                continue;
            }
            let distance_xy = distance_xy_sq.sqrt();
            let radial_speed = (tx * tvx + ty * tvy + tz * tvz) / distance;
            let longitude_speed = if distance_xy_sq > 0.0 {
                ((tx * tvy - ty * tvx) / distance_xy_sq).to_degrees()
            } else {
                0.0
            };
            let latitude_speed = if distance_xy > 0.0 {
                ((tvz - tz * radial_speed / distance) / distance_xy).to_degrees()
            } else {
                0.0
            };
            let latitude_ratio = (tz / distance).clamp(-1.0, 1.0);
            *point = [
                ty.atan2(tx).to_degrees().rem_euclid(360.0),
                latitude_ratio.asin().to_degrees(),
                distance,
                if want_speed { longitude_speed } else { 0.0 },
                if want_speed { latitude_speed } else { 0.0 },
                if want_speed { radial_speed } else { 0.0 },
            ];
        }
        Ok(out)
    }

    /// Node/apse 4-tuple at a UT epoch: converts
    /// through the bound Delta-T provider, then the TT entry above, so
    /// the native pair is evaluated at the converted epoch.
    pub fn nod_aps_ut(
        &self,
        body_id: i32,
        jd_ut: f64,
        method: i32,
        flags: i32,
    ) -> Result<[PositionResult; 4]> {
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        self.nod_aps_tt(body_id, jd_tt, method, flags)
    }

    /// Node/apse 4-tuple at a TT epoch for an explicit surface site
    /// : the geocentric set computed with `FLG_TOPOCTR`
    /// stripped (so the siteless missing-observer refusal never fires
    /// here), reduced by the topocentric parallax at the central epoch,
    /// then the entry equatorial rotation and the output-flag shaping
    ///
    /// `jd_ut1` is the explicit UT1 partner for the terrestrial-rotation
    /// series: the model reads it from its global timescale
    /// (`ts.tt_jd(jd_tt).ut1`), so callers pass the bundled default
    /// timescale value ([`crate::time::default_timescale_ut1_from_tt`]),
    /// never the bound model Delta-T chain (at 2024 the two differ by
    /// ~4.2 s, a ~2 km observer displacement). The
    /// caller-configured model TT and the independent frame/time TT
    /// stay distinct channels, like the projected zero-point clock.
    /// Attached native site acquisition obtains fresh scalar-Time UT1 itself
    /// at each central/apparent side and owns the associated Earth raw reads;
    /// the supplied partner remains the callback-absent stored channel's clock.
    /// `site` is `(lon_deg, lat_deg, height_m)`. Without
    /// `FLG_TOPOCTR` this delegates to [`nod_aps_tt`](Self::nod_aps_tt)
    /// and ignores the site. Projected sidereal selections compose
    /// with the topo arm, through the
    /// private `nod_aps_projected_topo_tt` composition entry.
    /// Output order, units, method flags and undefined-point sentinels follow
    /// [`Self::nod_aps_tt`]. SPEED uses its domain-aware geocentric rates;
    /// the site arm then carries them through the Cartesian velocity path at
    /// the central epoch. Site longitude is east-positive and latitude is
    /// north-positive. Errors: [`Error::NonFiniteArgument`] on a non-finite
    /// epoch or site component, plus the inner-call errors.
    pub fn nod_aps_topo_tt(
        &self,
        body_id: i32,
        jd_tt: f64,
        jd_ut1: f64,
        method: i32,
        flags: i32,
        site: (f64, f64, f64),
    ) -> Result<[PositionResult; 4]> {
        if (flags & constants::FLG_TOPOCTR) == 0 {
            return self.nod_aps_tt(body_id, jd_tt, method, flags);
        }
        if self.nodaps_projection_route(flags).is_some() {
            return self.nod_aps_projected_topo_tt(body_id, jd_tt, jd_ut1, method, flags, site);
        }
        self.nod_aps_plain_topo_tt(body_id, jd_tt, jd_ut1, method, flags, site)
    }

    /// Unprojected topocentric node/apse 4-tuple: the geocentric set
    /// computed with `FLG_TOPOCTR` stripped, reduced by the parallax
    /// at the central epoch, then the entry equatorial rotation and
    /// the output-flag shaping. Shared by [`nod_aps_topo_tt`](Self::nod_aps_topo_tt) and the
    /// `None`-matrix fallback of the projected-topo composition.
    fn nod_aps_plain_topo_tt(
        &self,
        body_id: i32,
        jd_tt: f64,
        jd_ut1: f64,
        method: i32,
        flags: i32,
        site: (f64, f64, f64),
    ) -> Result<[PositionResult; 4]> {
        let geo_flags = flags & !constants::FLG_TOPOCTR;
        let pts = self.nod_aps_inner_tt(body_id, jd_tt, method, geo_flags)?;
        let pts = self.apply_nodaps_topocentric(pts, jd_tt, jd_ut1, flags, site)?;
        let pts = self.apply_nodaps_equatorial(pts, jd_tt, flags)?;
        Ok(shape_nodaps_output(pts, flags))
    }

    /// Node/apse 4-tuple at a UT epoch for an explicit surface site
    /// : converts through the
    /// bound Delta-T provider (the matched-epoch isolation, like
    /// [`nod_aps_ut`](Self::nod_aps_ut)), then the TT topo entry above
    /// with the UT1 partner from the bundled default timescale at the
    /// converted epoch.
    pub fn nod_aps_topo_ut(
        &self,
        body_id: i32,
        jd_ut: f64,
        method: i32,
        flags: i32,
        site: (f64, f64, f64),
    ) -> Result<[PositionResult; 4]> {
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let jd_ut1 = crate::time::default_timescale_ut1_from_tt(jd_tt)?;
        self.nod_aps_topo_tt(body_id, jd_tt, jd_ut1, method, flags, site)
    }

    /// Rigid-rotation arm for a projected `FLG_SIDEREAL` request
    /// : the
    /// ecliptic set rotates into the target frame, or the equatorial
    /// set into the mean equator of t 0.
    fn nodaps_projection_route(&self, flags: i32) -> Option<NodapsProjectionRoute> {
        use crate::sidereal::SIDBIT_PROJECTION_SUPPRESS_MODES;
        if (flags & constants::FLG_SIDEREAL) == 0 {
            return None;
        }
        let selection = self.ayanamsha.as_ref()?.selection();
        let mode = strip_to_base_mode(selection.mode());
        let fixed = is_fixed_epoch_request(flags, mode);
        let bits = selection.bits();
        let sidbit = (bits & (constants::SIDBIT_ECL_T0 | constants::SIDBIT_SSY_PLANE)) != 0
            && !SIDBIT_PROJECTION_SUPPRESS_MODES.contains(&mode);
        if !fixed && !sidbit {
            return None;
        }
        if (flags & constants::FLG_EQUATORIAL) != 0 {
            // `SIDBIT_SSY_PLANE` alone leaves the equator untouched, like
            // the model: the scalar base arm serves it.
            if sidbit && !fixed && (bits & constants::SIDBIT_ECL_T0) == 0 {
                return None;
            }
            return Some(NodapsProjectionRoute::EquatorialT0);
        }
        Some(NodapsProjectionRoute::Ecliptic)
    }

    /// Projected node/apse 4-tuple at a TT epoch: the of-date set (or the
    /// J2000 set, when the caller set `FLG_J2000`) re-requested through
    /// this entry without `FLG_SIDEREAL` — so no wrapper recursion —
    /// with `FLG_NONUT` added and the representation bits stripped, then
    /// rotated rigidly into the target frame with the speeds. The
    /// rotation is defined on spherical degrees; shaping applies once
    /// below. Non-finite points pass through unrotated (the shared
    /// [`rotate_spherical_6tuple`](crate::sidereal::rotate_spherical_6tuple)
    /// NaN contract, same as the entry equatorial arm): the swept grids
    /// carry finite values only, and the interpolated-apsi sentinels stay
    /// outside the projected comparison scope with the lunar/minor
    /// points.
    ///
    /// Composition: without the caller's
    /// `FLG_J2000` the sub-request is the mean-of-date set and the
    /// rotation carries it into the target frame
    /// (`m_target @ ecliptic_of_t0(jd_tt).T`, equatorially
    /// `equatorial_epoch_matrix(t0) @ equatorial_epoch_matrix(jd_tt).T`);
    /// with it the J2000 reduction is composed into the rotation
    /// (`m_target`, equatorially `equatorial_epoch_matrix(t0)`). The
    /// `None` matrix arm falls back to the plain entry, like the
    /// the model `m_sidbit is None` arm.
    fn nod_aps_projected_tt(
        &self,
        body_id: i32,
        jd_tt: f64,
        method: i32,
        flags: i32,
    ) -> Result<[PositionResult; 4]> {
        use crate::sidereal::rotate_spherical_6tuple;
        let Some((sub_flags, rot)) = self.nodaps_projection_sub_and_rot(flags, jd_tt)? else {
            return self.nod_aps_plain_tt(body_id, jd_tt, method, flags);
        };
        let sub = self.nod_aps_tt(body_id, jd_tt, method, sub_flags)?;
        let mut out = sub;
        for pt in out.iter_mut() {
            *pt = rotate_spherical_6tuple(*pt, rot)?;
        }
        Ok(shape_nodaps_output(out, flags))
    }

    /// Projected topocentric node/apse 4-tuple at a TT epoch for an
    /// explicit surface site: the model
    /// re-requests the full entry — parallax included — through
    /// `_nodaps_sidereal_frame_projection`, so the sub-request carries
    /// `FLG_TOPOCTR` and the rigid rotation applies to the
    /// topocentric set. The geocentric set is computed with
    /// `FLG_TOPOCTR` stripped, reduced by the parallax at the central
    /// epoch under the sub-request flags (mean frame: `FLG_NONUT`
    /// added, `FLG_SIDEREAL` stripped, so the observer sits in the
    /// same frame the rotation expects), then the entry equatorial
    /// rotation (a no-op unless the equatorial projection applies),
    /// the rigid rotation, and the output-flag shaping once at the
    /// end. `jd_ut1`/`site` follow the
    /// [`nod_aps_topo_tt`](Self::nod_aps_topo_tt) convention. The
    /// `None`-matrix arm falls back to the plain topo entry, like the
    /// the model `m_sidbit is None` arm.
    fn nod_aps_projected_topo_tt(
        &self,
        body_id: i32,
        jd_tt: f64,
        jd_ut1: f64,
        method: i32,
        flags: i32,
        site: (f64, f64, f64),
    ) -> Result<[PositionResult; 4]> {
        use crate::sidereal::rotate_spherical_6tuple;
        let Some((sub_flags, rot)) = self.nodaps_projection_sub_and_rot(flags, jd_tt)? else {
            return self.nod_aps_plain_topo_tt(body_id, jd_tt, jd_ut1, method, flags, site);
        };
        let sub_flags = sub_flags | constants::FLG_TOPOCTR;
        let geo_sub = sub_flags & !constants::FLG_TOPOCTR;
        let pts = self.nod_aps_inner_tt(body_id, jd_tt, method, geo_sub)?;
        let pts = self.apply_nodaps_topocentric(pts, jd_tt, jd_ut1, sub_flags, site)?;
        let pts = self.apply_nodaps_equatorial(pts, jd_tt, sub_flags)?;
        let mut out = pts;
        for pt in out.iter_mut() {
            *pt = rotate_spherical_6tuple(*pt, rot)?;
        }
        Ok(shape_nodaps_output(out, flags))
    }

    /// Sub-request flags and rigid rotation for a projected
    /// `FLG_SIDEREAL` request: the sub-request is always spherical and
    /// non-sidereal — the rigid rotation is defined on that frame, and
    /// `FLG_XYZ` / `FLG_RADIANS` are representations of the rotated
    /// result — with `FLG_NONUT` added and `FLG_TOPOCTR` kept, so a
    /// topocentric caller carries the parallax inside the sub-request
    /// . Returns `None`
    /// when the SIDBIT matrix arm is missing, in which case the caller falls back to its plain entry.
    fn nodaps_projection_sub_and_rot(
        &self,
        flags: i32,
        jd_tt: f64,
    ) -> Result<Option<(i32, [[f64; 3]; 3])>> {
        use crate::positions::mat3_mul;
        use crate::sidereal::{
            J2000_TT, ecl_t0_epoch_jd, ecliptic_of_t0_matrix, equatorial_epoch_matrix,
            fixed_epoch_matrices, fixed_epoch_t0_jd, sidbit_ecliptic_matrix,
            ssy_plane_zero_point_deg,
        };
        use crate::time::default_timescale_tt_from_ut1;
        let route = self
            .nodaps_projection_route(flags)
            .unwrap_or(NodapsProjectionRoute::Ecliptic);
        let aya = self
            .ayanamsha
            .as_ref()
            .ok_or(Error::UnsupportedFlagCombination {
                entry: "nod_aps",
                flags,
            })?;
        let selection = aya.selection();
        let mode = strip_to_base_mode(selection.mode());
        let bits = selection.bits();
        let fixed = is_fixed_epoch_request(flags, mode);
        let transpose = |m: [[f64; 3]; 3]| {
            [
                [m[0][0], m[1][0], m[2][0]],
                [m[0][1], m[1][1], m[2][1]],
                [m[0][2], m[1][2], m[2][2]],
            ]
        };
        // The sub-request is always spherical and non-sidereal: the
        // rigid rotation is defined on that frame, and `FLG_XYZ` /
        // `FLG_RADIANS` are representations of the rotated result.
        let base = (flags
            & !(constants::FLG_SIDEREAL | constants::FLG_RADIANS | constants::FLG_XYZ))
            | constants::FLG_NONUT;
        let j2000 = (flags & constants::FLG_J2000) != 0;
        let (sub_flags, rot) = match route {
            NodapsProjectionRoute::Ecliptic => {
                let m_target = if fixed {
                    fixed_epoch_matrices(mode)?.1
                } else {
                    let (user_t0_tt, _) = aya.user_pair_tt()?;
                    let valens_tt = aya.valens_epoch_tt()?;
                    // Zero-point clock: the defining-epoch mean
                    // ayanamsha through the bundled default timescale,
                    // never the bound model Delta-T chain — at 1955 the
                    // two clocks differ by ~0.1 s, which the Method-B
                    // accumulation turns into a ~5 e-11 deg longitude
                    // shift (discriminating probe 2). Pinned by
                    // `sidbit_zero_point_independent_clock_matches_released_words`.
                    let zero_at = |t0_as_ut: f64| -> Result<f64> {
                        let t0_tt = default_timescale_tt_from_ut1(t0_as_ut)?;
                        Ok(aya.mean_tt_with(t0_tt, false, false)?.0)
                    };
                    let (t0_jd, zero_point) = if (bits & constants::SIDBIT_ECL_T0) != 0 {
                        let t0 = ecl_t0_epoch_jd(mode, user_t0_tt, valens_tt);
                        (t0, zero_at(t0)?)
                    } else {
                        (J2000_TT, ssy_plane_zero_point_deg(zero_at(J2000_TT)?)?)
                    };
                    match sidbit_ecliptic_matrix(bits, t0_jd, zero_point)? {
                        Some(m) => m,
                        None => {
                            return Ok(None);
                        }
                    }
                };
                if j2000 {
                    (base, m_target)
                } else {
                    let m = mat3_mul(m_target, transpose(ecliptic_of_t0_matrix(jd_tt)?))?;
                    (base & !constants::FLG_J2000, m)
                }
            }
            NodapsProjectionRoute::EquatorialT0 => {
                let t0_eq = if fixed {
                    fixed_epoch_t0_jd(mode).ok_or(Error::AyanamshaAnchorUnserved { mode })?
                } else {
                    let (user_t0_tt, _) = aya.user_pair_tt()?;
                    ecl_t0_epoch_jd(mode, user_t0_tt, aya.valens_epoch_tt()?)
                };
                if j2000 {
                    (base, equatorial_epoch_matrix(t0_eq)?)
                } else {
                    let m = mat3_mul(
                        equatorial_epoch_matrix(t0_eq)?,
                        transpose(equatorial_epoch_matrix(jd_tt)?),
                    )?;
                    (base & !constants::FLG_J2000, m)
                }
            }
        };
        Ok(Some((sub_flags, rot)))
    }

    /// Mean ayanamsha in degrees for a sidereal node/apse epoch
    /// : always the mean value. A native public feed reads the
    /// original observation Time's global UT1; callback-absent sealed callers
    /// retain their configured inverse. Dispatch warnings are dropped like the
    /// swept sidereal arms: the served selections are warning-free.
    fn nodaps_mean_ayanamsha_deg(&self, jd_tt: f64, flags: i32) -> Result<f64> {
        let Some(aya) = &self.ayanamsha else {
            return Err(Error::UnsupportedFlagCombination {
                entry: "nod_aps",
                flags,
            });
        };
        let ut1 = match self.nodaps_ut1_at {
            Some(at) => at(jd_tt)?,
            None => aya.ut1_for_tt(jd_tt)?,
        };
        let stripped = (flags & !constants::FLG_J2000) | constants::FLG_NONUT;
        self.nodaps_ayanamsha_at_ut1(aya, ut1, stripped)
    }

    /// True ayanamsha in degrees for the topocentric observer rotation
    /// . Native public acquisition uses a fresh scalar-TT Time inverse
    /// for every observer frame epoch, separately from the original UT Time;
    /// callback-absent sealed callers retain the configured inverse.
    fn nodaps_true_ayanamsha_deg(&self, jd_tt: f64, flags: i32) -> Result<f64> {
        let Some(aya) = &self.ayanamsha else {
            return Err(Error::UnsupportedFlagCombination {
                entry: "nod_aps",
                flags,
            });
        };
        let ut1 = if self.nodaps_ut1_at.is_some() {
            // Observer _to_result_frame creates a fresh scalar tt_jd Time,
            // even at the original UT entry's central scalar TT observable.
            crate::time::default_timescale_ut1_from_tt(jd_tt)?
        } else {
            aya.ut1_for_tt(jd_tt)?
        };
        let stripped = flags & !constants::FLG_J2000;
        self.nodaps_ayanamsha_at_ut1(aya, ut1, stripped)
    }

    /// Flag-aware scalar ayanamsha at the source Time's UT1 observable.
    /// Native observation TT comes from the independent global timescale;
    /// configured Delta-T remains attached to USER_UT/Valens defining epochs.
    /// Stored callers retain their configured observation conversion. Anchor
    /// aberration/deflection and mean/true toggles follow the supplied flags.
    fn nodaps_ayanamsha_at_ut1(
        &self,
        aya: &SealedAyanamshaServing<'_>,
        ut1: f64,
        flags: i32,
    ) -> Result<f64> {
        if self.nodaps_ut1_at.is_none() {
            return Ok(aya.for_flags_ut(ut1, flags)?.0);
        }
        let global_tt = crate::time::default_timescale_tt_from_ut1(ut1)?;
        let (noaberr, nogdefl, use_true) = crate::sidereal::resolve_ayanamsha_flags(flags);
        let value = if use_true {
            aya.true_tt_with(global_tt, noaberr, nogdefl)?
        } else {
            aya.mean_tt_with(global_tt, noaberr, nogdefl)?
        };
        Ok(value.0)
    }

    /// Analytic mean-lunar conic rates before translating or rotating its frame.
    /// The radius equation p/(1±e cos u) is differentiated by the chain rule;
    /// the common nutation cancels from u. Rates use the actual inner model
    /// clock and the selected nutation polynomial's local derivative.
    fn mean_moon_nodaps_rates(
        &self,
        jd_tt: f64,
        method: i32,
        flags: i32,
    ) -> Result<Option<[[f64; 3]; 4]>> {
        if !prefer_mean(method)
            || flags
                & (constants::FLG_HELCTR
                    | constants::FLG_BARYCTR
                    | constants::FLG_J2000
                    | constants::FLG_SIDEREAL)
                != 0
        {
            return Ok(None);
        }
        // An external native callback owns its own model clock; an actual
        // stored nutation channel identifies the bound LEB clock convention.
        if self.nodaps_calc_ut.is_some()
            && self
                .serving
                .provider()
                .nutation_smooth_interval_tt(jd_tt)?
                .is_none()
        {
            return Ok(None);
        }
        let ut1 = match self.nodaps_ut1_at {
            Some(clock) => clock(jd_tt)?,
            None => crate::time::default_timescale_ut1_from_tt(jd_tt)?,
        };
        let kernel_tt = ut1 + (self.deltat_at)(ut1)?;
        let shape_flags = (flags | constants::FLG_SPEED | constants::FLG_NONUT)
            & !(constants::FLG_RADIANS
                | constants::FLG_SPEED3
                | constants::FLG_XYZ
                | constants::FLG_EQUATORIAL
                | constants::FLG_TOPOCTR);
        let (node, apogee) = if let Some(calc_ut) = self.nodaps_calc_ut {
            (
                calc_ut(ut1, constants::MEAN_NODE, shape_flags)?.0,
                calc_ut(ut1, constants::MEAN_APOG, shape_flags)?.0,
            )
        } else {
            let node = crate::lunar::mean_lunar_node_state(
                kernel_tt,
                crate::lunar::LunarCoverage::Unenforced,
                false,
            )?;
            let apogee = crate::lunar::mean_lunar_apogee_state(
                kernel_tt,
                crate::lunar::LunarCoverage::Unenforced,
                false,
            )?;
            (
                [node.0, node.1, node.2, node.3, node.4, node.5],
                [apogee.0, apogee.1, apogee.2, apogee.3, apogee.4, apogee.5],
            )
        };
        let nutation_rate = if flags & constants::FLG_NONUT != 0 {
            0.0
        } else {
            match self.serving.provider().nutation_rate_tt(kernel_tt)? {
                Some((longitude, _)) => longitude.to_degrees(),
                None => {
                    let (start, end) = self
                        .serving
                        .provider()
                        .nutation_smooth_interval_tt(kernel_tt)?
                        .unwrap_or((f64::NEG_INFINITY, f64::INFINITY));
                    let (offsets, step) = nodaps_derivative_stencil(kernel_tt, start, end)?;
                    let nutation = self.serving.nutation_provider();
                    let central = nutation(kernel_tt)?.0;
                    let mut delta = [0.0; 4];
                    for (value, offset) in delta.iter_mut().zip(offsets) {
                        *value = nutation(kernel_tt + offset * step)?.0 - central;
                    }
                    nodaps_stencil_derivative(delta, offsets, step).to_degrees()
                }
            }
        };
        // Differentiate the small clock offsets. Subtracting full Julian days
        // would quantize this factor at tens of microseconds per sample.
        let clock_offset = |tt: f64| -> Result<f64> {
            let ut = crate::time::default_timescale_ut1_from_tt(tt)?;
            Ok((self.deltat_at)(ut)? - crate::smh2016::bundled_deltat_seconds(tt)? / 86400.0)
        };
        let h = NODAPS_SPEED_STEP_DAYS;
        let mut clock_delta = [0.0; 4];
        for (value, offset) in clock_delta.iter_mut().zip([-2.0, -1.0, 1.0, 2.0]) {
            *value = clock_offset(jd_tt + offset * h)?;
        }
        let clock_factor = 1.0 + nodaps_stencil_derivative(clock_delta, [-2.0, -1.0, 1.0, 2.0], h);
        let apog_radius = crate::lunar::MEAN_APOGEE_DISTANCE_AU;
        let peri_radius = 2.0 * MEAN_LUNAR_DISTANCE_AU - apog_radius;
        let a = (peri_radius + apog_radius) * 0.5;
        let e = (apog_radius - peri_radius) / (apog_radius + peri_radius);
        let p = a * (1.0 - e * e);
        let u = (node[0] - apogee[0] - 180.0).to_radians();
        let u_rate = (node[3] - apogee[3]).to_radians() * clock_factor;
        let radial_numerator = p * e * u.sin() * u_rate;
        let ascending_radius_rate = radial_numerator / (1.0 + e * u.cos()).powi(2);
        let descending_radius_rate = -radial_numerator / (1.0 - e * u.cos()).powi(2);
        let node_rate = (node[3] + nutation_rate) * clock_factor;
        let apogee_rate = (apogee[3] + nutation_rate) * clock_factor;
        Ok(Some([
            [node_rate, 0.0, ascending_radius_rate],
            [node_rate, 0.0, descending_radius_rate],
            [apogee_rate, -apogee[4] * clock_factor, 0.0],
            [apogee_rate, apogee[4] * clock_factor, 0.0],
        ]))
    }

    /// Unshaped spherical node/apse computation behind both entries.
    /// Speeds difference the unshaped triples; shaping happens once at the entry.
    fn nod_aps_inner_tt(
        &self,
        body_id: i32,
        jd_tt: f64,
        method: i32,
        flags: i32,
    ) -> Result<[PositionResult; 4]> {
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        if body_id == constants::INTP_APOG || body_id == constants::INTP_PERG {
            // Undefined interpolated distances are NaN, rather than the
            // zero-distance sentinels exempt from the observer requirement.
            // A centered request has already stripped TOPOCTR at the entry.
            if (flags & constants::FLG_TOPOCTR) != 0 {
                return Err(Error::ConfigurationMissing {
                    missing_config: "observer_location",
                    suggestion: "call set_topo(lon, lat, alt) first",
                });
            }
            return Ok(intp_undefined_points((flags & constants::FLG_SPEED) != 0));
        }
        if !is_nodaps_planet(body_id) && !is_minor_nodaps_body(body_id) {
            return Err(Error::InvalidOrbitalObject {
                entry: "nod_aps",
                detail: format!("nodes/apsides for body id {body_id} are not implemented."),
            });
        }
        if (flags & constants::FLG_TOPOCTR) != 0 {
            return Err(Error::UnsupportedFlagCombination {
                entry: "nod_aps",
                flags,
            });
        }
        // Plain `FLG_SIDEREAL` reduces to the scalar mean-ayanamsha arm
        // below (projected selections never reach the inner call: the
        // entry routes them through the rigid frame projection first);
        // a missing ayanamsha serving keeps the loud refusal, which is
        // also the production binding behavior.
        let want_sidereal = (flags & constants::FLG_SIDEREAL) != 0;
        if want_sidereal && self.ayanamsha.is_none() {
            return Err(Error::UnsupportedFlagCombination {
                entry: "nod_aps",
                flags,
            });
        }
        if (flags & constants::FLG_SPEED) != 0 {
            let base = flags & !constants::FLG_SPEED;
            if body_id == constants::MOON {
                if let Some(rates) = self.mean_moon_nodaps_rates(jd_tt, method, base)? {
                    let mut points = self.nod_aps_inner_tt(body_id, jd_tt, method, base)?;
                    for (point, rates) in points.iter_mut().zip(rates) {
                        point[3..].copy_from_slice(&rates);
                    }
                    return Ok(points);
                }
            }
            // The central evaluation identifies every polynomial support used
            // by the geometry, including retarded deflectors and lunar clocks.
            // Clear tracking before propagating an error or taking side samples.
            self.smooth_support.set(Some(NodapsSmoothSupport {
                epoch: jd_tt,
                start: f64::NEG_INFINITY,
                end: f64::INFINITY,
            }));
            let central = self.nod_aps_inner_tt(body_id, jd_tt, method, base);
            let support = self
                .smooth_support
                .take()
                .expect("central support is active");
            let now = central?;
            let (offsets, step) = nodaps_derivative_stencil(jd_tt, support.start, support.end)?;
            let mut samples = [now; 4];
            for (sample, offset) in samples.iter_mut().zip(offsets) {
                *sample = self.nod_aps_inner_tt(body_id, jd_tt + offset * step, method, base)?;
            }
            return Ok(std::array::from_fn(|point| {
                let position = now[point];
                let mut result = position;
                for component in 0..3 {
                    let delta = samples.map(|sample| {
                        let delta = sample[point][component] - position[component];
                        if component == 0 && delta > 180.0 {
                            delta - 360.0
                        } else if component == 0 && delta < -180.0 {
                            delta + 360.0
                        } else {
                            delta
                        }
                    });
                    result[component + 3] = nodaps_stencil_derivative(delta, offsets, step);
                }
                result
            }));
        }

        // Effective nutation for this computation: `FLG_NONUT` selects
        // the mean frame; otherwise the
        // native IAU 2006/2000 A pair at this epoch, never
        // the stored LEB series.
        let nut_eff = if (flags & constants::FLG_NONUT) != 0 {
            None
        } else {
            let (dpsi, deps) = crate::erfa::nut06a(jd_tt)?;
            Some((dpsi, deps))
        };
        let frame = NodapsFrame::for_date(jd_tt, nut_eff)?;
        if body_id == constants::MOON {
            return self.moon_nod_aps_tt(jd_tt, method, flags, nut_eff);
        }
        if is_minor_nodaps_body(body_id) {
            // Attached native acquisition uses the genuine public TT place;
            // callback-absent stored callers retain analytic Keplerian states.
            // Both reduce heliocentric ICRS through the shared conic tail.
            return self.minor_nod_aps_tt(body_id, jd_tt, method, flags, frame, nut_eff);
        }
        self.planet_nod_aps_tt(body_id, jd_tt, method, flags, frame, nut_eff)
    }

    /// Major-planet (and Sun-mirror) node/apse points at a TT epoch,
    /// unshaped.
    fn planet_nod_aps_tt(
        &self,
        body_id: i32,
        jd_tt: f64,
        method: i32,
        flags: i32,
        frame: NodapsFrame,
        nut_eff: Option<(f64, f64)>,
    ) -> Result<[PositionResult; 4]> {
        let request = NodapsPointRequest::new(body_id, jd_tt, method, flags);
        let mirror_sun = body_id == constants::SUN;
        let elem_body = if mirror_sun {
            constants::EARTH
        } else {
            body_id
        };
        // Resolve both observer aliases before evaluating the Sun. Range
        // errors here name the failing support body, unlike elements.
        let (sun_state, earth_state) = self.nodaps_observer_states(jd_tt)?;
        let (sun_p, _) = sun_state;
        let (earth_p, earth_v) = earth_state;
        // The public vector path evaluates the target even for a mean
        // fit. Reuse observer reads where the target is Sun or Earth.
        let target_state = match body_id {
            constants::SUN => sun_state,
            constants::EARTH => earth_state,
            _ => {
                if let Some(lookup) = self.nodaps_target_lookup {
                    lookup(body_id)?;
                }
                self.nodaps_icrs_state(body_id, jd_tt)?
            }
        };
        let earth_helio_ecl = frame.icrs_to_ecliptic([
            earth_p[0] - sun_p[0],
            earth_p[1] - sun_p[1],
            earth_p[2] - sun_p[2],
        ]);
        let earth_bary_ecl = frame.icrs_to_ecliptic(earth_p);
        let earth_vel_ecl = frame.icrs_to_ecliptic(earth_v);
        let observer = NodapsAcquiredObserver {
            frame,
            nut_eff,
            sun_p,
            earth_p,
            earth_helio_ecl,
            earth_bary_ecl,
            earth_vel_ecl,
        };

        let use_mean = prefer_mean(method) && !is_minor_nodaps_body(body_id);
        let mean_elements = if use_mean {
            mean_orbital_elements(elem_body, jd_tt)
        } else {
            None
        };
        // The orbit shape plus whether it refers to the mean equinox
        // (nutation-in-longitude alignment applies to mean shapes only).
        // Earth and the solar mirror fit the Earth-Moon barycenter,
        // while the observer remains the Earth's geocenter.
        let (shape, mean_frame) = match mean_elements {
            Some(el) => (
                mean_shape_from_elements(el[0], el[1], el[2], el[3], el[4]),
                true,
            ),
            None => {
                let (r_icrs, v_icrs) = if self.bar_mode(body_id, method) {
                    // SSB-relative vectors for the planets beyond
                    // Jupiter.
                    target_state
                } else {
                    // Heliocentric vectors.
                    let (tp, tv) = if elem_body == constants::EARTH {
                        if self.nodaps_target_lookup.is_some() {
                            self.emb_icrs_state(jd_tt)?
                        } else {
                            let moon = self.nodaps_icrs_state(constants::MOON, jd_tt)?;
                            emb_icrs_state_from_channels(earth_state, moon)
                        }
                    } else {
                        target_state
                    };
                    let (_, sun_v) = sun_state;
                    (
                        [tp[0] - sun_p[0], tp[1] - sun_p[1], tp[2] - sun_p[2]],
                        [tv[0] - sun_v[0], tv[1] - sun_v[1], tv[2] - sun_v[2]],
                    )
                };
                let gm = planet_system_gm(elem_body).unwrap_or(GM_SUN);
                (
                    osculating_shape_from_ecliptic(
                        frame.icrs_to_ecliptic(r_icrs),
                        frame.icrs_to_ecliptic(v_icrs),
                        gm,
                    ),
                    false,
                )
            }
        };
        // Mean elements describe a heliocentric ellipse even when the
        // method also carries OSCU_BAR. Only an actual osculating fit
        // can use barycentric vectors (including Pluto's mean fallback).
        let bar_mode = !mean_frame && self.bar_mode(body_id, method);
        self.finish_osculating_points(
            request,
            NodapsOrbitFit {
                shape,
                mean_frame,
                mirror_sun,
                bar_mode,
            },
            observer,
        )
    }

    /// Resolve Sun and geocenter Earth before either epoch read, then acquire
    /// them in that order. Native aliases use the injected resolver; sealed
    /// states retain the original channel/error policy. Returns raw SSB ICRS
    /// position/velocity pairs in AU and AU/day, without public calc guards.
    #[allow(clippy::type_complexity)]
    fn nodaps_observer_states(
        &self,
        jd_tt: f64,
    ) -> Result<(([f64; 3], [f64; 3]), ([f64; 3], [f64; 3]))> {
        for support in [constants::SUN, constants::EARTH] {
            if let Some(lookup) = self.nodaps_target_lookup {
                lookup(support)?;
            } else {
                crate::positions::serving_channel_gate(self.serving.provider(), support)
                    .map_err(|error| orbital_element_state_error(support, error))?;
            }
        }
        Ok((
            self.nodaps_icrs_state(constants::SUN, jd_tt)?,
            self.nodaps_icrs_state(constants::EARTH, jd_tt)?,
        ))
    }

    /// Checked ICRS state for a node/apse acquisition. Missing aliases
    /// and out-of-range reads retain the actual stored support identity;
    /// the public orbital-element requested-target remapping does not apply.
    fn nodaps_icrs_state(&self, body: i32, jd_tt: f64) -> Result<([f64; 3], [f64; 3])> {
        self.icrs_state(body, jd_tt)
            .map_err(|error| match self.nodaps_state_error {
                Some(at) => at(error),
                None => orbital_element_state_error(body, error),
            })
    }

    /// Compose the conic tail over already acquired observer operands.
    /// Request flags have one interpretation; the fit independently records
    /// its geometric origin. Deflector support callbacks are composed here,
    /// while point projection belongs to the borrowed projection context.
    /// The J2000 rule, scalar mean-ayanamsha read and zero-speed shaping
    /// follow projection in that order.
    fn finish_osculating_points(
        &self,
        request: NodapsPointRequest,
        fit: NodapsOrbitFit,
        observer: NodapsAcquiredObserver,
    ) -> Result<[PositionResult; 4]> {
        let NodapsPointRequest {
            body_id,
            jd_tt,
            method,
            flags,
            helio,
            bary,
            want_sidereal,
            want_equatorial,
        } = request;
        let NodapsOrbitFit {
            shape,
            mean_frame,
            mirror_sun,
            bar_mode,
        } = fit;
        let NodapsAcquiredObserver {
            frame,
            nut_eff,
            sun_p,
            earth_p,
            earth_helio_ecl,
            earth_bary_ecl,
            earth_vel_ecl,
        } = observer;
        let centered = helio || bary;
        let focal = (method & constants::NODBIT_FOPOINT) != 0;
        let asc = orbit_node_3d(shape, true);
        let dsc = orbit_node_3d(shape, false);
        let peri = orbit_point_3d(shape, 0.0);
        let aphe = if focal {
            focal_point_3d(shape)
        } else {
            orbit_point_3d(shape, std::f64::consts::PI)
        };

        let apply_aberr = !centered
            && !mirror_sun
            && (flags & (constants::FLG_TRUEPOS | constants::FLG_NOABERR)) == 0;
        let apply_defl = !centered
            && !mirror_sun
            && (flags & (constants::FLG_TRUEPOS | constants::FLG_NOGDEFL)) == 0;
        // Observer selection. OSCU_BAR
        // points are SSB-relative, so a barycentric request needs no
        // subtraction; otherwise the heliocentric Earth (or the
        // barycentric Sun/Earth for the bar/bary-centered arms).
        let observer_ecl = if bar_mode && bary {
            None
        } else if bar_mode && centered {
            Some(frame.icrs_to_ecliptic(sun_p))
        } else if bar_mode {
            Some(earth_bary_ecl)
        } else if bary {
            // FLG_BARYCTR on heliocentric orbit points re-references
            // each point to the SSB (the point gains the Sun's
            // barycentric offset), so the observer sits at -r_sun.
            Some(frame.icrs_to_ecliptic([-sun_p[0], -sun_p[1], -sun_p[2]]))
        } else if centered {
            None
        } else {
            Some(earth_helio_ecl)
        };
        let mean_align = if mean_frame {
            nut_eff.map(|(dpsi, _)| dpsi)
        } else {
            None
        };
        let deflect = |id: i32, jd: f64| {
            self.note_body_support(id, jd)?;
            self.serving
                .provider()
                .geometric_state(id, jd)
                .map_err(|error| match self.nodaps_state_error {
                    Some(at) => at(error),
                    None => orbital_element_state_error(id, error),
                })
                .map(|st| st.pos())
        };
        let mut projector = NodapsProjection::new(
            observer_ecl,
            mirror_sun && !centered,
            apply_aberr,
            apply_defl,
            earth_vel_ecl,
            earth_p,
            sun_p,
            jd_tt,
            frame,
            mean_align,
            &deflect,
        )?;
        let optional_deflect = |id: i32, jd: f64| {
            self.note_body_support(id, jd)?;
            self.nodaps_deflect_at
                .expect("optional deflector is attached")(id, jd)
        };
        if self.nodaps_deflect_at.is_some() {
            projector = projector.with_optional_deflector_lookup(&optional_deflect);
        }
        let mut pts =
            projector.project_orbit_points([asc.0, dsc.0, peri.0, aphe.0], body_id, mirror_sun)?;
        if (flags & constants::FLG_J2000) != 0 {
            pts = apply_nodaps_j2000_rule(pts, observer_ecl, mirror_sun, jd_tt)?;
        }
        if want_sidereal && !want_equatorial {
            // sidereal arm: subtract the mean
            // ayanamsha from the longitudes. `FLG_EQUATORIAL` suppresses
            // the subtraction entirely (SID+EQU == EQU alone; the
            // equatorial rotation runs at the entry). The zero-distance
            // sentinel slots (Sun/Earth nodes) are skipped, like the
            // guard. Per-epoch evaluation covers
            // the `FLG_SPEED` stencil sides, like the model recursion.
            let aya_deg = self.nodaps_mean_ayanamsha_deg(jd_tt, flags)?;
            for pt in pts.iter_mut() {
                if pt.2 != 0.0 {
                    pt.0 = (pt.0 - aya_deg).rem_euclid(360.0);
                }
            }
        }
        Ok([
            [pts[0].0, pts[0].1, pts[0].2, 0.0, 0.0, 0.0],
            [pts[1].0, pts[1].1, pts[1].2, 0.0, 0.0, 0.0],
            [pts[2].0, pts[2].1, pts[2].2, 0.0, 0.0, 0.0],
            [pts[3].0, pts[3].1, pts[3].2, 0.0, 0.0, 0.0],
        ])
    }

    /// Whether the osculating fit runs barycentrically: `NODBIT_OSCU_BAR`
    /// for Saturn and outward.
    fn bar_mode(&self, body_id: i32, method: i32) -> bool {
        (method & constants::NODBIT_OSCU_BAR) != 0
            && matches!(
                body_id,
                constants::SATURN | constants::URANUS | constants::NEPTUNE | constants::PLUTO
            )
    }

    /// Lunar-theory node/apogee 6-tuple for the Moon branch
    /// : the mean
    /// points from the analytic Delaunay theory, the true node from
    /// the stored LEB channel, the osculating apogee from the sealed
    /// LEB Moon/Earth states through the true-ecliptic-of-date
    /// geometry.
    /// `jd_kernel` is the caller-chosen theory epoch (the lunar clock
    /// for fast-path bodies, TT for the slow-path apogee). `source_ut1`
    /// is the independent global-timescale request epoch used by the
    /// true-node coverage policy and slow fallback; `flags`
    /// are the stripped `calc_flags` (`FLG_NONUT` retained, so the
    /// of-date longitude is mean-of-date when requested).
    fn nodaps_lunar_point(
        &self,
        body_id: i32,
        jd_kernel: f64,
        source_ut1: f64,
        flags: i32,
    ) -> Result<PositionResult> {
        self.note_body_support(body_id, jd_kernel)?;
        self.note_nutation_support(jd_kernel)?;
        let nutation = self.serving.nutation_provider();
        let dpsi_rad_at = &|jd: f64| nutation(jd).map(|(dpsi, _)| dpsi);
        let nut_rate_at = &|jd: f64, half: f64| {
            if !half.is_finite() || half <= 0.0 {
                return Err(Error::NonFiniteArgument { value: half });
            }
            Ok((dpsi_rad_at(jd + half)? - dpsi_rad_at(jd - half)?).to_degrees() / (2.0 * half))
        };
        if body_id == constants::MEAN_NODE || body_id == constants::MEAN_APOG {
            crate::lunar::mean_lunar_state_with_providers(
                body_id,
                jd_kernel,
                flags,
                dpsi_rad_at,
                nut_rate_at,
            )
        } else if body_id == constants::TRUE_NODE {
            // Fast-path stored channel: the channel stores the
            // true-ecliptic longitude/latitude (the distance proxy is
            // replaced by the conic downstream, so only lon/lat travel
            // on). Under a nutation-free frame the stored nutation in
            // longitude is stripped with the same stored `dpsi` the
            // frame data serves. A missing channel falls back to the
            // live two-body geometry (the eclipse-layer precedent for
            // tiers without the record).
            self.nodaps_true_node_channel(jd_kernel, source_ut1, flags, dpsi_rad_at)
        } else {
            // Slow-path osculating apogee: live two-body
            // geometry in the equinox-based frame (above), with the native
            // ERFA nutation stripped under a nutation-free frame
            // . Speeds are zero (the caller strips `FLG_SPEED`;
            // rates come from the outer stencil).
            let (r, v) = self.nodaps_true_ecliptic_state(jd_kernel)?;
            let (lon, lat, dist) = crate::lunar::true_lilith_from_state(r, v)?;
            let lon = if crate::lunar::no_nutation_frame(flags) {
                let (dpsi, _) = crate::erfa::nut06a(jd_kernel)?;
                (lon - dpsi.to_degrees()).rem_euclid(360.0)
            } else {
                lon
            };
            Ok([lon, lat, dist, 0.0, 0.0, 0.0])
        }
    }

    /// True-ecliptic-of-date Moon-minus-Earth state at fresh scalar TT.
    ///
    /// An attached native feed resolves Earth before Moon and evaluates the
    /// reversed negative Earth constituents before the forward Moon chain,
    /// retaining fresh scalar-TT parts and per-segment AU/AU-day conversion.
    /// Without that feed the stored provider evaluates the relative polynomial
    /// directly when its records align. Split provider domains and missing/range/frame
    /// errors remain distinct from public calculation coverage. The analytic
    /// IAU 2006/2000 A frame combines P 03 precession at TDB, IAU 2000A nutation,
    /// bias and true obliquity. Returns ecliptic AU and AU/day for the conic
    /// reducer. Its native ERFA nutation removal is a separate operation.
    fn nodaps_true_ecliptic_state(&self, jd_tt: f64) -> Result<crate::lunar::EclipticMoonState> {
        self.note_body_support(constants::EARTH, jd_tt)?;
        self.note_body_support(constants::MOON, jd_tt)?;
        use crate::positions::reference_frames;
        if let Some(state_at) = self.nodaps_moon_state_at {
            return state_at(jd_tt);
        }
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        if let Some(lookup) = self.nodaps_target_lookup {
            // The scalar-TT perigee constructor resolves Earth before Moon,
            // while the relative vector still reads Moon before Earth.
            lookup(constants::EARTH)?;
            lookup(constants::MOON)?;
        }
        let (relative, velocity) =
            self.relative_icrs_state(constants::MOON, constants::EARTH, jd_tt)?;
        let tdb = jd_tt + reference_frames::tdb_minus_tt_seconds(jd_tt, 0.0)? / 86400.0;
        let eps = (reference_frames::mean_obliquity_arcsec(tdb)? / 3600.0).to_radians();
        let (psi, deps) = reference_frames::iau2000a_tenth_microarcsec(jd_tt)?;
        let scale = reference_frames::FRAME_TENTH_USEC_2_RAD;
        let nut = reference_frames::nutation_matrix(eps, eps + deps * scale, psi * scale)?;
        let equinox =
            reference_frames::equinox_matrix(nut, reference_frames::precession_matrix(tdb)?)?;
        let frame = reference_frames::ecliptic_date_matrix(equinox, eps, deps * scale)?;
        let r = reference_frames::matvec3(frame, relative)?;
        let v = reference_frames::matvec3(frame, velocity)?;
        Ok(((r[0], r[1], r[2]), (v[0], v[1], v[2])))
    }

    /// Stored true-node channel read for the Moon branch: stored ecliptic-direct lon/lat
    /// at the kernel epoch with the stored nutation stripped under a
    /// nutation-free frame. Speeds are zero (the caller strips
    /// `FLG_SPEED`; rates come from the outer stencil). A missing
    /// channel can fall back to live two-body geometry over the sealed
    /// Moon/Earth states. A stored-channel miss first checks the active
    /// true-node coverage at `source_ut1`, the request UT1 JD; an excluded
    /// request returns [`Error::EphemerisRange`]. An admitted request
    /// whose fast-model epoch misses the channel can still fall through.
    /// Live geometry and the nutation strip use TT reconstructed from
    /// `source_ut1` by the global default timescale, independently of the
    /// caller-model fast epoch `jd_kernel`. Other provider errors propagate.
    fn nodaps_true_node_channel(
        &self,
        jd_kernel: f64,
        source_ut1: f64,
        flags: i32,
        dpsi_rad_at: &dyn Fn(f64) -> Result<f64>,
    ) -> Result<PositionResult> {
        let live = || -> Result<PositionResult> {
            if let Some((start, end)) = self.serving.provider().body_coverage(constants::TRUE_NODE)
            {
                if !(start <= source_ut1 && source_ut1 <= end) {
                    return Err(Error::EphemerisRange {
                        requested_jd: source_ut1,
                        start_jd: start,
                        end_jd: end,
                        body: Some(constants::TRUE_NODE),
                    });
                }
            }
            let jd_slow = crate::time::default_timescale_tt_from_ut1(source_ut1)?;
            let (r, v) = self.nodaps_true_ecliptic_state(jd_slow)?;
            let (lon, lat, dist) = crate::lunar::true_lunar_node_from_state(r, v)?;
            let lon = if crate::lunar::no_nutation_frame(flags) {
                let (dpsi, _) = crate::erfa::nut06a(jd_slow)?;
                (lon - dpsi.to_degrees()).rem_euclid(360.0)
            } else {
                lon
            };
            Ok([lon, lat, dist, 0.0, 0.0, 0.0])
        };
        let channel = match self
            .serving
            .provider()
            .geometric_state(constants::TRUE_NODE, jd_kernel)
        {
            Ok(channel) => channel,
            Err(Error::LebBodyNotFound { .. }) | Err(Error::LebOutOfRange { .. }) => {
                return live();
            }
            Err(other) => return Err(other),
        };
        // The base-tier record is `COORD_ECLIPTIC`-direct; any other
        // non-ICRS ecliptic-direct layout is equally readable, while an
        // ICRS or unknown layout means the record is not the
        // fast-path channel this read mirrors.
        if !matches!(
            channel.frame(),
            crate::positions::StateFrame::EclipticDirect
                | crate::positions::StateFrame::GeocentricEcliptic
        ) {
            return Err(Error::SealedChannelNotServed {
                body_id: constants::TRUE_NODE,
                coord_type: self
                    .serving
                    .provider()
                    .serving_coord_type(constants::TRUE_NODE)
                    .unwrap_or(0),
            });
        }
        let triple = channel.pos();
        if !(triple[0].is_finite() && triple[1].is_finite() && triple[2].is_finite()) {
            return Err(Error::NonFiniteArgument { value: triple[0] });
        }
        let lon = if crate::lunar::no_nutation_frame(flags) {
            (triple[0] - dpsi_rad_at(jd_kernel)?.to_degrees()).rem_euclid(360.0)
        } else {
            triple[0]
        };
        Ok([lon, triple[1], triple[2], 0.0, 0.0, 0.0])
    }

    /// Lunar node/apse points: node/apogee reads
    /// through the lunar-theory dispatch above (mean or osculating
    /// lunar-theory bodies by the model choice), perigee from the
    /// osculating ellipse or the mean identity. The composition after
    /// the reads is the shared tested core.
    fn moon_nod_aps_tt(
        &self,
        jd_tt: f64,
        method: i32,
        flags: i32,
        nut_eff: Option<(f64, f64)>,
    ) -> Result<[PositionResult; 4]> {
        let helio = (flags & constants::FLG_HELCTR) != 0;
        let bary = (flags & constants::FLG_BARYCTR) != 0 && !helio;
        let want_sidereal = (flags & constants::FLG_SIDEREAL) != 0;
        let want_equatorial = (flags & constants::FLG_EQUATORIAL) != 0;
        let mean = prefer_mean(method);
        let native_support = if let Some(lookup) = self.nodaps_target_lookup {
            let supports = self.nodaps_observer_states(jd_tt)?;
            match lookup(constants::MOON) {
                Ok(()) => {}
                // Moon has no source barycenter fallback; its missing target
                // alias returns four zero sentinels after observer acquisition.
                Err(Error::SpkUnregisteredBody {
                    body: constants::MOON,
                }) => {
                    return Ok([[0.0; 6]; 4]);
                }
                Err(other) => return Err(other),
            }
            self.nodaps_icrs_state(constants::MOON, jd_tt)?;
            Some(supports)
        } else if self.nodaps_calc_ut.is_some() {
            // Source-facing stored composition consumes all three geometric
            // supports before either public model point. The Moon sample is
            // mandatory even when mean elements later replace its shape.
            let supports = self.nodaps_observer_states(jd_tt)?;
            self.nodaps_icrs_state(constants::MOON, jd_tt)?;
            Some(supports)
        } else {
            None
        };
        let (node_body, apog_body) = if mean {
            (constants::MEAN_NODE, constants::MEAN_APOG)
        } else {
            (constants::TRUE_NODE, constants::OSCU_APOG)
        };
        // Spherical of-date places (output-format bits stripped, like
        // `calc_flags` mask; FLG_NONUT retained so the
        // of-date longitude is mean-of-date when requested).
        let calc_flags = flags
            & !constants::FLG_SPEED
            & !constants::FLG_RADIANS
            & !constants::FLG_XYZ
            & !constants::FLG_J2000
            & !constants::FLG_EQUATORIAL
            & !constants::FLG_SIDEREAL
            & !constants::FLG_HELCTR
            & !constants::FLG_BARYCTR
            & !constants::FLG_TOPOCTR;
        // Native public reads preserve calc_ut(Time.ut 1) as a distinct route.
        // The callback-absent sealed clock below follows fast_calc_ut:
        // its calc_ut reads inside the Moon branch evaluate the
        // theory at UT1 plus the caller Delta-T chain, where UT1 is the
        // bundled default timescale value, never the TT argument itself.
        // Retaining the independent UT1 conversion matters because the
        // caller Delta-T model can differ from the global timescale.
        let (node_pos, apog_pos) = if let Some(calc_ut) = self.nodaps_calc_ut {
            let jd_ut1 = self.nodaps_ut1_at.expect("native nodaps clock is attached")(jd_tt)?;
            // Native/SPK callbacks own their clocks and have no piecewise
            // LEB support. Obtain a model clock only while differentiating a
            // stored channel or stored nutation, never for a plain position.
            if self.smooth_support.get().is_some()
                && (self
                    .serving
                    .provider()
                    .smooth_interval_tt(node_body, jd_tt)?
                    .is_some()
                    || self
                        .serving
                        .provider()
                        .smooth_interval_tt(apog_body, jd_tt)?
                        .is_some()
                    || self
                        .serving
                        .provider()
                        .nutation_smooth_interval_tt(jd_tt)?
                        .is_some())
            {
                let jd_kernel = jd_ut1 + (self.deltat_at)(jd_ut1)?;
                self.note_body_support(node_body, jd_kernel)?;
                self.note_body_support(apog_body, jd_kernel)?;
                self.note_nutation_support(jd_kernel)?;
            }
            let node = calc_ut(jd_ut1, node_body, calc_flags)?.0;
            let apog = calc_ut(jd_ut1, apog_body, calc_flags)?.0;
            (node, apog)
        } else {
            let jd_ut1 = crate::time::default_timescale_ut1_from_tt(jd_tt)?;
            let dt_days = (self.deltat_at)(jd_ut1)?;
            if !dt_days.is_finite() {
                return Err(Error::NonFiniteArgument { value: dt_days });
            }
            let jd_eval = jd_ut1 + dt_days;
            let node_pos = self.nodaps_lunar_point(node_body, jd_eval, jd_ut1, calc_flags)?;
            // Backend split inside the model Moon branch: the node read
            // and the mean apogee read run through the sealed `calc_ut`
            // fast path (evaluation at the lunar clock above), while the
            // osculating apogee misses its LEB channel and falls through
            // to the slow `_calc_body` path at the TT epoch.
            let apog_epoch = if mean { jd_eval } else { jd_tt };
            let apog_pos = self.nodaps_lunar_point(apog_body, apog_epoch, jd_ut1, calc_flags)?;
            (node_pos, apog_pos)
        };
        // The mean lunar theory defines a fixed apogee radius. Preserve that
        // scalar through direction/frame rotations instead of redifferentiating
        // a rounded Cartesian norm. Translated helio/barycentric points still
        // acquire their actual, time-dependent distance in the composition.
        let apog_dist = if mean {
            crate::lunar::MEAN_APOGEE_DISTANCE_AU
        } else {
            apog_pos[2]
        };
        let peri_dist = if mean {
            2.0 * MEAN_LUNAR_DISTANCE_AU - apog_dist
        } else {
            // Osculating perigee distance from the instantaneous
            // ellipse at the TT epoch.
            let (r, v) = self.nodaps_true_ecliptic_state(jd_tt)?;
            crate::lunar::osculating_perigee_from_state(r, v)?.2
        };
        let ((sun_p, _), (earth_p, _)) = match native_support {
            Some(supports) => supports,
            None => (
                self.icrs_state(constants::SUN, jd_tt)?,
                self.icrs_state(constants::EARTH, jd_tt)?,
            ),
        };
        // Earth support vectors in the working of-date ecliptic (the
        // closure over the same frame).
        let frame = NodapsFrame::for_date(jd_tt, nut_eff)?;
        let earth_center_ecl = frame.icrs_to_ecliptic([
            earth_p[0] - sun_p[0],
            earth_p[1] - sun_p[1],
            earth_p[2] - sun_p[2],
        ]);
        let earth_bary_ecl = frame.icrs_to_ecliptic(earth_p);
        // Sidereal arm: the same mean ayanamsha
        // as the planet path, suppressed under `FLG_EQUATORIAL`. The
        // lunar points never hit the zero-distance sentinel, so no skip
        // applies here.
        let ayanamsha_deg = if want_sidereal && !want_equatorial {
            Some(self.nodaps_mean_ayanamsha_deg(jd_tt, flags)?)
        } else {
            None
        };
        moon_nodaps_points(
            node_pos[0],
            node_pos[1],
            apog_pos[0],
            apog_pos[1],
            apog_dist,
            peri_dist,
            (method & constants::NODBIT_FOPOINT) != 0,
            helio,
            bary,
            earth_center_ecl,
            earth_bary_ecl,
            (flags & constants::FLG_J2000) != 0,
            jd_tt,
            ayanamsha_deg,
        )
    }

    /// Unshaped minor-body node/apse points at scalar TT.
    ///
    /// Attached native acquisition reads Sun/Earth raw SPK support, then the
    /// genuine configured public TT calculation with HELCTR|J2000|TRUEPOS|SPEED.
    /// Its degree/AU/per-day spherical state is converted through the J2000
    /// frame transpose into heliocentric ICRS AU/AU-day. Callback-absent stored
    /// callers use the analytic Keplerian state. Public model/provider clocks,
    /// configuration, output guards and typed errors belong to that callback.
    /// Every method uses osculating elements for minors; OSCU_BAR adds the SSB
    /// Sun offset for orbits beyond TRANS_JOVIAN_A_AU before the shared conic
    /// reduction. Missing registrations/models and range faults propagate;
    /// no zero state or alternate provider is substituted.
    fn minor_nod_aps_tt(
        &self,
        body_id: i32,
        jd_tt: f64,
        method: i32,
        flags: i32,
        frame: NodapsFrame,
        nut_eff: Option<(f64, f64)>,
    ) -> Result<[PositionResult; 4]> {
        let request = NodapsPointRequest::new(body_id, jd_tt, method, flags);
        let native_support = if self.nodaps_target_lookup.is_some() {
            Some(self.nodaps_observer_states(jd_tt)?)
        } else {
            None
        };
        let (r_icrs, v_icrs) = if native_support.is_some() {
            let calc_flags = constants::FLG_HELCTR
                | constants::FLG_J2000
                | constants::FLG_TRUEPOS
                | constants::FLG_SPEED;
            let position = self.model_calc_tt(jd_tt, body_id, calc_flags)?.0;
            minor_icrs_state_from_j2000_position(position)?
        } else {
            minor_helio_icrs_state(body_id, jd_tt)?
        };
        // Trans-jovian re-reference under `NODBIT_OSCU_BAR`: the heliocentric state gains the Sun's
        // barycentric offset. The classifier runs on the heliocentric
        // state with plain solar GM (no mass-ratio entry exists for
        // minor bodies).
        let mut bar_mode = false;
        let (r_fit, v_fit) = if (method & constants::NODBIT_OSCU_BAR) != 0
            && osculating_semi_major_axis(
                r_icrs[0], r_icrs[1], r_icrs[2], v_icrs[0], v_icrs[1], v_icrs[2], GM_SUN,
            ) > TRANS_JOVIAN_A_AU
        {
            bar_mode = true;
            let (sp, sv) = match native_support {
                Some((sun, _)) => sun,
                None => self.icrs_state(constants::SUN, jd_tt)?,
            };
            (
                [r_icrs[0] + sp[0], r_icrs[1] + sp[1], r_icrs[2] + sp[2]],
                [v_icrs[0] + sv[0], v_icrs[1] + sv[1], v_icrs[2] + sv[2]],
            )
        } else {
            (r_icrs, v_icrs)
        };
        let shape = osculating_shape_from_ecliptic(
            frame.icrs_to_ecliptic(r_fit),
            frame.icrs_to_ecliptic(v_fit),
            GM_SUN,
        );
        // Support states (SSB-relative ICRS, single fetch each).
        let ((sun_p, _), (earth_p, earth_v)) = match native_support {
            Some(supports) => supports,
            None => (
                self.icrs_state(constants::SUN, jd_tt)?,
                self.icrs_state(constants::EARTH, jd_tt)?,
            ),
        };
        let earth_helio_ecl = frame.icrs_to_ecliptic([
            earth_p[0] - sun_p[0],
            earth_p[1] - sun_p[1],
            earth_p[2] - sun_p[2],
        ]);
        let earth_bary_ecl = frame.icrs_to_ecliptic(earth_p);
        let earth_vel_ecl = frame.icrs_to_ecliptic(earth_v);
        let observer = NodapsAcquiredObserver {
            frame,
            nut_eff,
            sun_p,
            earth_p,
            earth_helio_ecl,
            earth_bary_ecl,
            earth_vel_ecl,
        };
        self.finish_osculating_points(
            request,
            NodapsOrbitFit {
                shape,
                mean_frame: false,
                mirror_sun: false,
                bar_mode,
            },
            observer,
        )
    }
}

/// Heliocentric geometric ICRS state (AU, AU/day) for a minor body
/// : the native epoch-blended
/// Keplerian propagation shaped to TRUE-of-date heliocentric spherical,
/// stripped to mean-of-date under the nutation-free J2000 request,
/// rated by the model 1-second stencil, precessed to the mean J2000
/// ecliptic with the state precession map, shaped to Cartesian with the
/// analytic derivative, and rotated to ICRS by the transpose of the
/// J2000 frame matrix — in the same ICRS frame as the planetary states used by the shared reduction. Unknown ids are
/// [`Error::UnknownBody`]; degenerate distances are
/// [`Error::NonFiniteArgument`], and every other provider miss
/// propagates fail-closed.
fn minor_helio_icrs_state(body_id: i32, jd_tt: f64) -> Result<([f64; 3], [f64; 3])> {
    // Mean-of-date heliocentric spherical at one epoch through the
    // native Keplerian chain.
    let mean_of_date_at = |jd: f64| -> Result<(f64, f64, f64)> {
        let (dpsi, _) = crate::erfa::nut06a(jd)?;
        let (lon_true, lat_true, dist) =
            crate::minor_bodies::calc_minor_body_heliocentric(body_id, jd, false, None, dpsi)?;
        Ok((
            (lon_true - dpsi.to_degrees()).rem_euclid(360.0),
            lat_true,
            dist,
        ))
    };
    let (lon, lat, dist) = mean_of_date_at(jd_tt)?;
    // the model 1-second central stencil with the longitude
    // wrap fold.
    const STENCIL_DT_DAYS: f64 = 1.0 / 86400.0;
    let (lon_prev, lat_prev, dist_prev) = mean_of_date_at(jd_tt - STENCIL_DT_DAYS)?;
    let (lon_next, lat_next, dist_next) = mean_of_date_at(jd_tt + STENCIL_DT_DAYS)?;
    let mut dlon_wrapped = lon_next - lon_prev;
    if dlon_wrapped > 180.0 {
        dlon_wrapped -= 360.0;
    } else if dlon_wrapped < -180.0 {
        dlon_wrapped += 360.0;
    }
    let dlon = dlon_wrapped / (2.0 * STENCIL_DT_DAYS);
    let dlat = (lat_next - lat_prev) / (2.0 * STENCIL_DT_DAYS);
    let ddist = (dist_next - dist_prev) / (2.0 * STENCIL_DT_DAYS);
    // Mean-J2000-ecliptic state.
    let (lon_j2000, lat_j2000, dlon_j2000, dlat_j2000) =
        crate::hypothetical::precess_ecliptic_state(lon, lat, dlon, dlat, jd_tt, true)?;
    minor_icrs_state_from_j2000_position([
        lon_j2000, lat_j2000, dist, dlon_j2000, dlat_j2000, ddist,
    ])
}

/// Convert a public mean-J2000 ecliptic spherical state into heliocentric
/// ICRS AU/AU-day. The six values are degrees, AU and per-day rates. A finite
/// positive distance is required; the checked Cartesian conversion rejects
/// non-finite state components. The transpose uses the same fixed frame and
/// written association as the source public minor-state helper.
fn minor_icrs_state_from_j2000_position(position: PositionResult) -> Result<([f64; 3], [f64; 3])> {
    let dist = position[2];
    if !dist.is_finite() || dist <= 0.0 {
        return Err(Error::NonFiniteArgument { value: dist });
    }
    let cart = ecliptic_spherical_to_cartesian_state(position);
    for value in cart {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    // Mean-J2000-ecliptic -> ICRS is the transpose of the ICRS ->
    // ecliptic frame matrix.
    let m = icrs_to_ecl_j2000_matrix();
    let to_icrs = |ax: f64, ay: f64, az: f64| -> [f64; 3] {
        [
            m[0][0] * ax + m[1][0] * ay + m[2][0] * az,
            m[0][1] * ax + m[1][1] * ay + m[2][1] * az,
            m[0][2] * ax + m[1][2] * ay + m[2][2] * az,
        ]
    };
    Ok((
        to_icrs(cart[0], cart[1], cart[2]),
        to_icrs(cart[3], cart[4], cart[5]),
    ))
}

/// True for bodies whose elements reduce from the position pipeline's
/// geometric place: curated minor
/// bodies, numbered asteroids (`AST_OFFSET <= ipl < FIXSTAR_OFFSET`,
/// inclusive lower end so the pipeline's typed `UnknownBodyError` for
/// "asteroid 0" propagates), and fictitious ids.
fn is_unserved_model_body(body_id: i32) -> bool {
    matches!(
        body_id,
        constants::CHIRON
            | constants::PHOLUS
            | constants::CERES
            | constants::PALLAS
            | constants::JUNO
            | constants::VESTA
    ) || (constants::AST_OFFSET..constants::FIXSTAR_OFFSET).contains(&body_id)
        || (constants::CUPIDO..=constants::WALDEMATH).contains(&body_id)
}

/// Samples per eccentric-anomaly axis on the coarse distance-extremum grid.
pub const ORBIT_EXTREMUM_COARSE_N: usize = 1440;

/// Shrinking-window iterations for each local distance-extremum refinement.
pub const ORBIT_EXTREMUM_REFINE_ITERS: usize = 25;

/// Equally spaced samples per anomaly axis in one local refinement window,
/// including both endpoints.
pub const ORBIT_EXTREMUM_REFINE_WIDTH: usize = 9;

/// Heliocentric osculating ellipse behind the extremum search.
/// Fields are private with accessors under the construction invariant.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OrbitEllipse {
    a_au: f64,
    e: f64,
    inc_rad: f64,
    node_rad: f64,
    arg_peri_rad: f64,
}

/// Earth/Moon mass ratio behind the sealed EMB composition.
pub const EMB_EARTH_MOON_MASS_RATIO: f64 = 81.3005691;

/// EMB composition weight `1 / (Earth/Moon mass ratio + 1)`.
pub fn emb_composition_weight() -> f64 {
    1.0 / (EMB_EARTH_MOON_MASS_RATIO + 1.0)
}

/// True for the geocentric fictitious constructions: White Moon Selena and the Waldemath dark
/// moon, whose osculating ellipse is taken about the Earth.
pub fn is_fict_geocentric_body(body_id: i32) -> bool {
    body_id == constants::WHITE_MOON || body_id == constants::WALDEMATH
}

/// Compose an Earth-Moon-barycenter state from acquired Earth and Moon states.
///
/// Both inputs are `(position, velocity)` in AU and AU/day, relative to the
/// same barycentric origin and ICRS axes at the same epoch. This pure weighted
/// sum performs no acquisition or validation; nonfinite arithmetic propagates.
/// The caller owns all channel errors and epoch/frame consistency checks.
pub fn emb_icrs_state_from_channels(
    earth: ([f64; 3], [f64; 3]),
    moon: ([f64; 3], [f64; 3]),
) -> ([f64; 3], [f64; 3]) {
    let w = emb_composition_weight();
    let ((ep, ev), (mp, mv)) = (earth, moon);
    (
        [
            ep[0] + (mp[0] - ep[0]) * w,
            ep[1] + (mp[1] - ep[1]) * w,
            ep[2] + (mp[2] - ep[2]) * w,
        ],
        [
            ev[0] + (mv[0] - ev[0]) * w,
            ev[1] + (mv[1] - ev[1]) * w,
            ev[2] + (mv[2] - ev[2]) * w,
        ],
    )
}

/// Compatibility no-op: `FLG_HELCTR` is ignored for the Sun and the
/// Moon (`hel == geo`; head).
pub fn strip_helio_for_sun_moon(body_id: i32, flags: i32) -> i32 {
    if body_id == constants::SUN || body_id == constants::MOON {
        flags & !constants::FLG_HELCTR
    } else {
        flags
    }
}

/// `FLG_BARYCTR` interior normalization: an explicit `FLG_BARYCTR`
/// request for a body whose heliocentric `a` does not exceed
/// [`TRANS_JOVIAN_A_AU`] degrades to `FLG_HELCTR` so the true-distance
/// slot and the element branch stay mutually heliocentric.
pub fn normalize_barycentric_flag(body_id: i32, flags: i32, a_helio_au: f64) -> i32 {
    // the model predicate `not (a_helio > 6.0)` (NaN normalizes,
    // +inf keeps the request): spelled without a negated comparison so
    // the NaN arm stays explicit under clippy.
    let interior = a_helio_au <= TRANS_JOVIAN_A_AU || a_helio_au.is_nan();
    if (flags & constants::FLG_BARYCTR) != 0
        && body_id != constants::SUN
        && body_id != constants::MOON
        && interior
    {
        (flags & !constants::FLG_BARYCTR) | constants::FLG_HELCTR
    } else {
        flags
    }
}

impl<'a, P> SealedOrbitalServing<'a, P>
where
    P: StateProvider,
{
    /// Maximum, minimum and current geometric position component
    /// . Returns `(max, min, current)`; the extrema are distances in AU.
    /// `current` preserves public `calc(TT, body, flags | TRUEPOS)[0][2]`
    /// . With spherical output it is distance in AU; under XYZ it is
    /// the signed z component in the requested ecliptic/equatorial frame.
    /// Combined public representation flags retain their shaping semantics,
    /// so this component must not always be interpreted as a distance.
    ///
    /// Served arms: heliocentric/barycentric extremes from the sealed
    /// elements (Sun: the model zero tuple; Earth: the EMB ellipse),
    /// the geocentric Moon from its lunar elements, the geocentric
    /// fictitious points from their geocentric pipeline place, the
    /// geocentric Sun/Earth/planet extremes from the sealed EMB
    /// composition, and the geometric true-distance slot.
    /// Minor/fictitious positions use the attached public calculation feed;
    /// unknown models and missing support propagate typed errors. Without a
    /// feed they retain the stored-channel serving policy. Geocentric two-orbit
    /// extrema use the bounded approximation in [`two_ellipse_min_max`], not
    /// simultaneous positions on an N-body trajectory or a certified global
    /// optimum. Unbound/unsupported ellipse arms can return zero extrema while
    /// retaining the separately calculated current component.
    pub fn orbit_max_min_true_distance_tt(
        &self,
        body_id: i32,
        jd_tt: f64,
        flags: i32,
    ) -> Result<(f64, f64, f64)> {
        validate_orbital_object("orbit_max_min_true_distance", body_id, true)?;
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let mut eff = strip_helio_for_sun_moon(body_id, flags);
        if (eff & constants::FLG_BARYCTR) != 0
            && body_id != constants::SUN
            && body_id != constants::MOON
        {
            // A typed source error is the caller's answer, not a reason
            // to silently downgrade an explicit FLG_BARYCTR request
            // ; genuinely unresolvable geometry falls back to
            // the heliocentric convention.
            match self.helio_semi_major_au(body_id, jd_tt) {
                Ok(a) => eff = normalize_barycentric_flag(body_id, eff, a),
                Err(
                    e @ (Error::EphemerisRange { .. }
                    | Error::UnknownBody { .. }
                    | Error::SealedChannelNotServed { .. }),
                ) => {
                    return Err(e);
                }
                Err(_) => eff = normalize_barycentric_flag(body_id, eff, 0.0),
            }
        }
        let true_dist = self.orbit_true_distance_tt(body_id, jd_tt, eff)?;
        if (eff & (constants::FLG_HELCTR | constants::FLG_BARYCTR)) != 0 {
            // The Sun has no elements: the zero
            // triple. Earth's heliocentric orbit is the EMB orbit
            // reduced with the Earth-system GM.
            if body_id == constants::SUN {
                return Ok((0.0, 0.0, true_dist));
            }
            let el = if body_id == constants::EARTH {
                self.earth_helio_elements_tt(jd_tt)?
            } else {
                self.elements_inner(body_id, jd_tt, eff)?
            };
            return Ok(match helio_distance_extremes(&el) {
                Some((max, min)) => (max, min, true_dist),
                None => (0.0, 0.0, true_dist),
            });
        }
        if body_id == constants::MOON {
            let el = self.moon_elements_tt(jd_tt)?;
            return Ok(match helio_distance_extremes(&el) {
                Some((max, min)) => (max, min, true_dist),
                None => (0.0, 0.0, true_dist),
            });
        }
        if is_fict_geocentric_body(body_id) {
            let el = self.calc_state_geocentric_elements(body_id, jd_tt, eff)?;
            return Ok(match helio_distance_extremes(&el) {
                Some((max, min)) => (max, min, true_dist),
                None => (0.0, 0.0, true_dist),
            });
        }
        // Sun/Earth acquire their EMB ellipse directly. Other targets first
        // acquire and validate their own ellipse, so a target error or an
        // unbound trajectory precedes unnecessary EMB support reads.
        if body_id == constants::SUN {
            let emb = self.emb_elements_tt(jd_tt)?;
            let (max, min) = sun_extremes_from_emb_elements(&emb);
            return Ok((max, min, true_dist));
        }
        if body_id == constants::EARTH {
            let emb = self.emb_elements_tt(jd_tt)?;
            let (max, min) = earth_extremes_from_emb_elements(&emb);
            return Ok((max, min, true_dist));
        }
        let el = self.elements_inner(body_id, jd_tt, eff)?;
        if el[0] <= 0.0 || el[1] >= 1.0 {
            return Ok((0.0, 0.0, true_dist));
        }
        let emb = self.emb_elements_tt(jd_tt)?;
        Ok(match geocentric_two_ellipse_extremes(&el, &emb)? {
            Some((max, min)) => (max, min, true_dist),
            None => (0.0, 0.0, true_dist),
        })
    }

    /// Maximum/minimum distances and current position component at a UT1 epoch.
    ///
    /// Convert through the bound model Delta-T callback, then use
    /// [`Self::orbit_max_min_true_distance_tt`]. Its result order, AU extrema,
    /// bounded-ellipse limitations and zero sentinels apply unchanged. The third
    /// value is a distance only for spherical output; XYZ retains signed z in
    /// the requested frame. Object validation precedes clock conversion.
    pub fn orbit_max_min_true_distance_ut(
        &self,
        body_id: i32,
        jd_ut: f64,
        flags: i32,
    ) -> Result<(f64, f64, f64)> {
        validate_orbital_object("orbit_max_min_true_distance", body_id, true)?;
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        self.orbit_max_min_true_distance_tt(body_id, jd_tt, flags)
    }

    /// Heliocentric semi-major axis behind the `FLG_BARYCTR` interior
    /// normalization: Earth's heliocentric orbit is the EMB orbit
    /// reduced with the Earth-system GM; every other id reduces through `elements_inner`.
    fn helio_semi_major_au(&self, body_id: i32, jd_tt: f64) -> Result<f64> {
        if body_id == constants::EARTH {
            Ok(self.earth_helio_elements_tt(jd_tt)?[0])
        } else {
            Ok(self.elements_inner(body_id, jd_tt, constants::FLG_HELCTR)?[0])
        }
    }

    /// Osculating heliocentric ecliptic elements of the Earth.
    /// The same reduction is exposed by the public Earth element entry
    /// when no Astronomical Almanac cumulative-mass flag is selected.
    fn earth_helio_elements_tt(&self, jd_tt: f64) -> Result<[f64; 50]> {
        let ((bp, bv), (sp, sv)) = if let Some(target_lookup) = self.orbital_target_lookup {
            // The source Earth element branch resolves EMB before the Sun,
            // then evaluates the Sun before the EMB. Keep legacy reads below.
            target_lookup(constants::EARTH)?;
            target_lookup(constants::SUN)?;
            let sun = self.icrs_state(constants::SUN, jd_tt)?;
            let emb = self.emb_icrs_state(jd_tt)?;
            (emb, sun)
        } else {
            let emb = self.emb_icrs_state(jd_tt)?;
            let sun = self.icrs_state(constants::SUN, jd_tt)?;
            (emb, sun)
        };
        let p = [bp[0] - sp[0], bp[1] - sp[1], bp[2] - sp[2]];
        let v = [bv[0] - sv[0], bv[1] - sv[1], bv[2] - sv[2]];
        let m = icrs_to_ecl_j2000_matrix();
        let pos = mat3_vec3(m, p)?;
        let vel = mat3_vec3(m, v)?;
        let gm = planet_system_gm(constants::EARTH).unwrap_or(GM_SUN);
        Ok(orbital_elements_from_ecliptic_state(
            [pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]],
            gm,
            constants::EARTH,
            jd_tt,
            Some(GM_SUN),
        ))
    }

    /// EMB ICRS state (AU, AU/day): an attached direct provider feed takes
    /// precedence; otherwise the sealed Earth/Moon composition implements
    /// `leb_vector._eval_earth_barycenter`. Neither path substitutes the
    /// geocenter for the Earth-Moon system barycenter.
    fn emb_icrs_state(&self, jd_tt: f64) -> Result<([f64; 3], [f64; 3])> {
        self.note_body_support(constants::EARTH, jd_tt)?;
        self.note_body_support(constants::MOON, jd_tt)?;
        if let Some(emb_at) = self.orbital_emb_at {
            let state = emb_at(jd_tt)?;
            return self.checked_orbital_icrs_state(constants::EARTH, jd_tt, state);
        }
        let earth = self.icrs_state(constants::EARTH, jd_tt)?;
        let moon = self.icrs_state(constants::MOON, jd_tt)?;
        Ok(emb_icrs_state_from_channels(earth, moon))
    }

    /// Osculating heliocentric ecliptic elements of the EMB.
    fn emb_elements_tt(&self, jd_tt: f64) -> Result<[f64; 50]> {
        let ((bp, bv), (sp, sv)) = if let Some(target_lookup) = self.orbital_target_lookup {
            // _emb_osculating_ecliptic_elements resolves both aliases before
            // either epoch read: Sun, EMB, then Sun.at(Time), EMB.at(Time).
            target_lookup(constants::SUN)?;
            target_lookup(constants::EARTH)?;
            let sun = self.icrs_state(constants::SUN, jd_tt)?;
            let emb = self.emb_icrs_state(jd_tt)?;
            (emb, sun)
        } else {
            let emb = self.emb_icrs_state(jd_tt)?;
            let sun = self.icrs_state(constants::SUN, jd_tt)?;
            (emb, sun)
        };
        let p = [bp[0] - sp[0], bp[1] - sp[1], bp[2] - sp[2]];
        let v = [bv[0] - sv[0], bv[1] - sv[1], bv[2] - sv[2]];
        let m = icrs_to_ecl_j2000_matrix();
        let pos = mat3_vec3(m, p)?;
        let vel = mat3_vec3(m, v)?;
        Ok(orbital_elements_from_ecliptic_state(
            [pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]],
            GM_SUN,
            constants::EARTH,
            jd_tt,
            None,
        ))
    }

    /// Elements from the pipeline geometric *geocentric* place. Caller flags do not change the center, frame or mass:
    /// ORBEL_AA is a heliocentric convention and BARYCTR is not forwarded.
    /// Public model acquisition errors propagate without a substituted state.
    fn calc_state_geocentric_elements(
        &self,
        body_id: i32,
        jd_tt: f64,
        flags: i32,
    ) -> Result<[f64; 50]> {
        let geo_flags = constants::FLG_J2000 | constants::FLG_TRUEPOS | constants::FLG_SPEED;
        let (pos, _) = self.model_calc_tt(jd_tt, body_id, geo_flags)?;
        if !pos[2].is_finite() || pos[2] <= 0.0 {
            return Err(Error::NonFiniteArgument { value: pos[2] });
        }
        let cart = ecliptic_spherical_to_cartesian_state(pos);
        let gm_earth = GM_SUN / SUN_EARTH_MASS_RATIO;
        let _ = flags;
        Ok(orbital_elements_from_ecliptic_state(
            [cart[0], cart[1], cart[2], cart[3], cart[4], cart[5]],
            gm_earth,
            body_id,
            jd_tt,
            None,
        ))
    }

    /// Geometric true-distance slot: the `pos[2]` value, so output
    /// shaping quirks are preserved. Typed provider misses propagate
    /// fail-closed.
    fn orbit_true_distance_tt(&self, body_id: i32, jd_tt: f64, flags: i32) -> Result<f64> {
        if jd_tt <= 0.0 {
            return Ok(0.0);
        }
        let calc_flags = flags | constants::FLG_TRUEPOS;
        let (pos, _) = if let Some(calc_tt) = self.current_calc_tt {
            calc_tt(jd_tt, body_id, calc_flags)?
        } else if is_unserved_model_body(body_id) {
            self.model_calc_tt(jd_tt, body_id, calc_flags)?
        } else {
            serving_compute_tt(self.serving, jd_tt, body_id, calc_flags)?
        };
        if !pos[2].is_finite() {
            return Err(Error::NonFiniteArgument { value: pos[2] });
        }
        Ok(pos[2])
    }
}

#[cfg(test)]
mod model_element_tests {
    use super::*;
    use std::cell::RefCell;

    /// Each model read must use the attached public calculation feed. This
    /// provider makes an accidental stored target/support acquisition visible.
    struct CallbackOnlyProvider;
    impl StateProvider for CallbackOnlyProvider {
        fn provider_name(&self) -> &'static str {
            "test/public-model-callback"
        }
        fn has_body(&self, _: i32) -> bool {
            panic!("model callback must not read a stored target or center")
        }
        fn body_coverage(&self, _: i32) -> Option<(f64, f64)> {
            panic!("model callback owns its coverage")
        }
        fn serving_coord_type(&self, _: i32) -> Option<u32> {
            panic!("model callback owns its frame")
        }
        fn eval_raw(&self, _: i32, _: f64) -> Result<(Vec<f64>, Vec<f64>)> {
            panic!("model callback owns its state")
        }
        fn eval_raw_two_part(
            &self,
            _: i32,
            _: crate::types::TwoPartEpoch,
        ) -> Result<(Vec<f64>, Vec<f64>)> {
            panic!("model callback owns its state")
        }
        fn serving_entry_snapshot(&self, _: i32) -> Result<crate::positions::EntrySnapshot> {
            panic!("model callback owns its entry")
        }
    }

    #[test]
    fn all_target_current_feed_preserves_ordinary_public_shaping_and_priority() {
        let provider = CallbackOnlyProvider;
        let nut =
            |_: f64| -> Result<(f64, f64)> { panic!("current callback owns frame acquisition") };
        let serving = SealedServing::new(&provider, &nut);
        let delta = |_: f64| -> Result<f64> { panic!("TT current does not convert its epoch") };
        let calls = RefCell::new(Vec::new());
        let calc = |tt, body, flags| {
            calls.borrow_mut().push((tt, body, flags));
            Ok(([1.0, 2.0, -0.003, 0.0, 0.0, 0.0], flags))
        };
        let bound = SealedOrbitalServing::new(&serving, &delta).with_current_calc_tt(&calc);
        let flags = constants::FLG_XYZ | constants::FLG_EQUATORIAL;
        let current = bound
            .orbit_true_distance_tt(constants::MARS, 2451545.125, flags)
            .unwrap();
        assert_eq!(current.to_bits(), (-0.003_f64).to_bits());
        assert_eq!(
            *calls.borrow(),
            vec![(2451545.125, constants::MARS, flags | constants::FLG_TRUEPOS)]
        );

        let refused = |_: f64, _: i32, _: i32| -> Result<(PositionResult, i32)> {
            Err(Error::ConfigurationMissing {
                missing_config: "observer_location",
                suggestion: "test public refusal",
            })
        };
        let bound = bound.with_current_calc_tt(&refused);
        assert!(matches!(
            bound.orbit_max_min_true_distance_tt(
                constants::MARS,
                2451545.125,
                constants::FLG_TOPOCTR
            ),
            Err(Error::ConfigurationMissing {
                missing_config: "observer_location",
                ..
            })
        ));
        // The panic provider proves that a current-position error precedes
        // any target, Sun or EMB state evaluation in the extrema route.
    }

    #[test]
    fn orbital_target_alias_feed_refuses_before_center_geometry() {
        let provider = CallbackOnlyProvider;
        let nut = |_: f64| -> Result<(f64, f64)> { panic!("no stored frame read") };
        let serving = SealedServing::new(&provider, &nut);
        let delta = |_: f64| -> Result<f64> { panic!("no TT conversion") };
        let raw = |_: i32, _: f64| -> Result<crate::positions::GeometricState> {
            panic!("target alias must precede center evaluation")
        };
        let emb = |_: f64| -> Result<crate::positions::GeometricState> {
            panic!("no EMB read after target alias failure")
        };
        let lookup = |body| Err(Error::SpkUnregisteredBody { body });
        let bound =
            SealedOrbitalServing::new(&serving, &delta).with_orbital_geometry(&raw, &emb, &lookup);
        assert!(matches!(
            bound.orbital_elements_tt(constants::MARS, 2451545.125, 0),
            Err(Error::SpkUnregisteredBody {
                body: constants::MARS
            })
        ));
        assert!(matches!(
            bound.orbital_elements_tt(constants::SUN, 2451545.125, 0),
            Err(Error::InvalidOrbitalObject { .. })
        ));
    }

    #[test]
    fn native_emb_alias_order_precedes_any_epoch_acquisition() {
        let provider = CallbackOnlyProvider;
        let nut = |_: f64| -> Result<(f64, f64)> { panic!("no stored frame read") };
        let serving = SealedServing::new(&provider, &nut);
        let delta = |_: f64| -> Result<f64> { panic!("no TT conversion") };
        let raw = |_: i32, _: f64| -> Result<crate::positions::GeometricState> {
            panic!("both aliases must precede an epoch read")
        };
        let emb = |_: f64| -> Result<crate::positions::GeometricState> {
            panic!("missing EMB alias must refuse before geometry")
        };
        let calls = RefCell::new(Vec::new());
        let lookup = |body| {
            calls.borrow_mut().push(body);
            if body == constants::EARTH {
                Err(Error::SpkBodyNotFound { target: 3 })
            } else {
                Ok(())
            }
        };
        let bound =
            SealedOrbitalServing::new(&serving, &delta).with_orbital_geometry(&raw, &emb, &lookup);
        assert!(matches!(
            bound.emb_elements_tt(2451545.125),
            Err(Error::SpkBodyNotFound { target: 3 })
        ));
        assert_eq!(*calls.borrow(), vec![constants::SUN, constants::EARTH]);
    }

    /// An exact circular construction supplies an independent center-mass
    /// check: its radius is the fitted semi-major axis for the selected GM.
    fn circular_position(gm: f64, radius: f64) -> PositionResult {
        [
            0.0,
            0.0,
            radius,
            (gm / (radius * radius * radius)).sqrt().to_degrees(),
            0.0,
            0.0,
        ]
    }

    #[test]
    fn fictitious_reduction_uses_public_model_center_frame_and_mass() {
        let provider = CallbackOnlyProvider;
        let no_nutation = |_: f64| -> Result<(f64, f64)> { panic!("no stored frame read") };
        let serving = SealedServing::new(&provider, &no_nutation);
        let no_delta = |_: f64| -> Result<f64> { panic!("TT does not convert its epoch") };
        let calls = RefCell::new(Vec::new());
        let calc = |tt, body, flags| {
            calls.borrow_mut().push((tt, body, flags));
            if body == constants::NIBIRU {
                return Err(Error::UnknownBody { body_id: body });
            }
            let (gm, radius) = if is_fict_geocentric_body(body) {
                (moon_central_gm(), 0.002)
            } else {
                (GM_SUN, 40.0)
            };
            Ok((circular_position(gm, radius), flags))
        };
        let bound = SealedOrbitalServing::new(&serving, &no_delta).with_calc_tt(&calc);
        let tt = 2460400.0;
        let geometry_flags = constants::FLG_J2000 | constants::FLG_TRUEPOS | constants::FLG_SPEED;
        for body in [constants::WHITE_MOON, constants::WALDEMATH] {
            let plain = bound.orbital_elements_tt(body, tt, 0).unwrap();
            assert!((plain[0] - 0.002).abs() < 1e-15);
            assert!(plain[1] < 1e-12);
            let ignored = constants::FLG_BARYCTR
                | constants::FLG_ORBEL_AA
                | constants::FLG_XYZ
                | constants::FLG_RADIANS
                | constants::FLG_SIDEREAL;
            let repeated = bound.orbital_elements_tt(body, tt, ignored).unwrap();
            assert_eq!(plain.map(f64::to_bits), repeated.map(f64::to_bits));
            assert_eq!(*calls.borrow().last().unwrap(), (tt, body, geometry_flags));
        }
        for body in [constants::CUPIDO, constants::PROSERPINA] {
            let plain = bound.orbital_elements_tt(body, tt, 0).unwrap();
            assert!((plain[0] - 40.0).abs() < 1e-12);
            let bary = bound
                .orbital_elements_tt(body, tt, constants::FLG_BARYCTR)
                .unwrap();
            assert_eq!(plain.map(f64::to_bits), bary.map(f64::to_bits));
            let requested = constants::FLG_HELCTR
                | constants::FLG_TRUEPOS
                | constants::FLG_SPEED
                | if body == constants::PROSERPINA {
                    constants::FLG_NONUT
                } else {
                    constants::FLG_J2000
                };
            assert_eq!(*calls.borrow().last().unwrap(), (tt, body, requested));
            let aa = bound
                .orbital_elements_tt(body, tt, constants::FLG_ORBEL_AA)
                .unwrap();
            assert!(aa[0] < plain[0]);
            let solar_motion = (GM_SUN / aa[0].powf(3.0)).sqrt().to_degrees();
            assert_eq!(aa[11].to_bits(), solar_motion.to_bits());
        }
        assert!(matches!(
            bound.orbital_elements_tt(constants::NIBIRU, tt, 0),
            Err(Error::UnknownBody {
                body_id: constants::NIBIRU
            })
        ));
        let legacy = SealedOrbitalServing::new(&serving, &no_delta);
        assert!(matches!(
            legacy.orbital_elements_tt(constants::WHITE_MOON, tt, 0),
            Err(Error::UnknownBody {
                body_id: constants::WHITE_MOON
            })
        ));
    }

    #[test]
    fn model_extrema_reads_public_current_component_before_geocentric_reduction() {
        let provider = CallbackOnlyProvider;
        let no_nutation = |_: f64| -> Result<(f64, f64)> { panic!("no stored frame read") };
        let serving = SealedServing::new(&provider, &no_nutation);
        let no_delta = |_: f64| -> Result<f64> { panic!("TT does not convert its epoch") };
        let calls = RefCell::new(Vec::new());
        let calc = |tt, body, flags| {
            calls.borrow_mut().push((tt, body, flags));
            let state = if flags & constants::FLG_XYZ != 0 {
                [0.001, 0.002, -0.003, 0.0, 0.0, 0.0]
            } else {
                circular_position(moon_central_gm(), 0.002)
            };
            Ok((state, flags))
        };
        let bound = SealedOrbitalServing::new(&serving, &no_delta).with_calc_tt(&calc);
        let (max, min, current) = bound
            .orbit_max_min_true_distance_tt(constants::WHITE_MOON, 2460400.0, constants::FLG_XYZ)
            .unwrap();
        assert!((max - 0.002).abs() < 1e-15 && (min - 0.002).abs() < 1e-15);
        assert_eq!(current.to_bits(), (-0.003_f64).to_bits());
        assert_eq!(
            *calls.borrow(),
            vec![
                (
                    2460400.0,
                    constants::WHITE_MOON,
                    constants::FLG_XYZ | constants::FLG_TRUEPOS
                ),
                (
                    2460400.0,
                    constants::WHITE_MOON,
                    constants::FLG_J2000 | constants::FLG_TRUEPOS | constants::FLG_SPEED
                ),
            ]
        );
    }

    #[test]
    fn model_extrema_target_failure_and_unbound_orbit_precede_emb_reads() {
        let provider = CallbackOnlyProvider;
        let no_nutation = |_: f64| -> Result<(f64, f64)> { panic!("no stored frame read") };
        let serving = SealedServing::new(&provider, &no_nutation);
        let no_delta = |_: f64| -> Result<f64> { panic!("TT does not convert its epoch") };
        let unknown = |_, body, _| Err(Error::UnknownBody { body_id: body });
        let bound = SealedOrbitalServing::new(&serving, &no_delta).with_calc_tt(&unknown);
        assert!(matches!(
            bound.orbit_max_min_true_distance_tt(constants::CUPIDO, 0.0, 0),
            Err(Error::UnknownBody {
                body_id: constants::CUPIDO
            })
        ));
        // Speed twice the circular speed is above escape speed. The source
        // skips current-distance acquisition at TT <= 0, then returns its
        // undefined-ellipse zeros before reading Earth/Moon/Sun support.
        let unbound = |_, _, flags| {
            Ok((
                [0.0, 0.0, 1.0, (2.0 * GAUSS_K).to_degrees(), 0.0, 0.0],
                flags,
            ))
        };
        let bound = SealedOrbitalServing::new(&serving, &no_delta).with_calc_tt(&unbound);
        let result = bound
            .orbit_max_min_true_distance_tt(constants::CUPIDO, 0.0, 0)
            .unwrap();
        assert_eq!(result.0.to_bits(), 0.0_f64.to_bits());
        assert_eq!(result.1.to_bits(), 0.0_f64.to_bits());
        assert_eq!(result.2.to_bits(), 0.0_f64.to_bits());
    }

    #[test]
    fn orbital_timescale_injection_does_not_change_extension_or_legacy_clocks() {
        let provider = CallbackOnlyProvider;
        let no_nutation = |_: f64| -> Result<(f64, f64)> { panic!("no stored frame read") };
        let serving = SealedServing::new(&provider, &no_nutation);
        let delta = |_: f64| Ok(0.125);
        let calls = RefCell::new(Vec::new());
        let calc = |tt, _, flags| {
            calls.borrow_mut().push(tt);
            Ok((circular_position(moon_central_gm(), 0.002), flags))
        };
        let clock_calls = RefCell::new(Vec::new());
        let to_tt = |ut| {
            clock_calls.borrow_mut().push(ut);
            Ok(ut + 0.25)
        };
        let ut = 2451545.0;
        let legacy = SealedOrbitalServing::new(&serving, &delta).with_calc_tt(&calc);
        legacy
            .orbital_elements_ut(constants::WHITE_MOON, ut, 0)
            .unwrap();
        assert_eq!(*calls.borrow(), vec![ut + 0.125]);
        calls.borrow_mut().clear();
        let bound = legacy.with_orbital_ut1_to_tt(&to_tt);
        bound
            .orbital_elements_ut(constants::WHITE_MOON, ut, 0)
            .unwrap();
        assert_eq!(*calls.borrow(), vec![ut + 0.25]);
        assert_eq!(bound.jd_tt_for_ut(ut).unwrap(), ut + 0.125);
        calls.borrow_mut().clear();
        bound
            .orbit_max_min_true_distance_ut(constants::WHITE_MOON, ut, 0)
            .unwrap();
        assert_eq!(*calls.borrow(), vec![ut + 0.125, ut + 0.125]);
        assert!(matches!(
            bound.orbital_elements_ut(constants::SUN, ut, 0),
            Err(Error::InvalidOrbitalObject {
                entry: "get_orbital_elements_ut",
                ..
            })
        ));
        assert!(matches!(
            bound.orbital_elements_ut(constants::WHITE_MOON, f64::NAN, 0),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert_eq!(*clock_calls.borrow(), vec![ut]);
        let bad_clock = |_: f64| Ok(f64::INFINITY);
        let bad = bound.with_orbital_ut1_to_tt(&bad_clock);
        assert!(matches!(
            bad.orbital_elements_ut(constants::WHITE_MOON, ut, 0),
            Err(Error::NonFiniteArgument { .. })
        ));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn true_node_fallback_keeps_request_domain_and_global_slow_clock() {
        use crate::positions::SingleLeb;
        use crate::types::TwoPartEpoch;
        use sha2::{Digest, Sha256};

        struct NodeWindow<'a> {
            base: SingleLeb<'a>,
            window: Option<(f64, f64)>,
        }
        impl StateProvider for NodeWindow<'_> {
            fn provider_name(&self) -> &'static str {
                "test/partial-node"
            }
            fn has_body(&self, body: i32) -> bool {
                if body == constants::TRUE_NODE {
                    self.window.is_some()
                } else {
                    self.base.has_body(body)
                }
            }
            fn body_coverage(&self, body: i32) -> Option<(f64, f64)> {
                if body == constants::TRUE_NODE {
                    self.window
                } else {
                    self.base.body_coverage(body)
                }
            }
            fn serving_coord_type(&self, body: i32) -> Option<u32> {
                if self.has_body(body) {
                    self.base.serving_coord_type(body)
                } else {
                    None
                }
            }
            fn eval_raw(&self, body: i32, jd: f64) -> Result<(Vec<f64>, Vec<f64>)> {
                if body == constants::TRUE_NODE {
                    match self.window {
                        None => return Err(Error::LebBodyNotFound { body_id: body }),
                        Some((start, end)) if !(start <= jd && jd <= end) => {
                            return Err(Error::LebOutOfRange {
                                body_id: body,
                                jd,
                                start,
                                end,
                            });
                        }
                        _ => {}
                    }
                }
                self.base.eval_raw(body, jd)
            }
            fn eval_raw_two_part(
                &self,
                body: i32,
                epoch: TwoPartEpoch,
            ) -> Result<(Vec<f64>, Vec<f64>)> {
                self.eval_raw(body, epoch.value())
            }
            fn serving_entry_snapshot(&self, body: i32) -> Result<crate::positions::EntrySnapshot> {
                self.base.serving_entry_snapshot(body)
            }
        }
        // The scalar reference below was captured on the historical data-v3 distribution.
        // Do not apply that golden value to a different current-data image.
        let bytes = std::fs::read(crate::test_data::historical_core_path())
            .expect("the historical data-v3 regression fixture must exist");
        let hash: String = Sha256::digest(&bytes)
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect();
        assert_eq!(
            hash,
            "5d708bdbe3e799e0802ba575984e57a3c5e44720dbfa1b4a01cf826640e0cb82"
        );
        let index = crate::leb::FileIndex::decode(&bytes).unwrap();
        let jd_tt = 2460400.0;
        let request_ut1 = crate::time::default_timescale_ut1_from_tt(jd_tt).unwrap();
        let cfg = crate::time::DeltaTConfig {
            model: crate::time::DeltaTModel::EspenakMeeus,
            ..Default::default()
        };
        let delta_t = |ut| crate::time::deltat_with_config(ut, &cfg);
        let fast_tt = request_ut1 + delta_t(request_ut1).unwrap();
        let end = jd_tt - 0.0006;
        assert!(request_ut1 < end && end < fast_tt);
        for window in [None, Some((2396758.5, end)), Some((2396758.5, jd_tt - 1.0))] {
            let provider = NodeWindow {
                base: SingleLeb::new(&index, &bytes),
                window,
            };
            let nutation = |jd| index.eval_nutation(&bytes, jd);
            let serving = SealedServing::new(&provider, &nutation);
            let orbital = SealedOrbitalServing::new(&serving, &delta_t);
            let dpsi = |jd| nutation(jd).map(|(psi, _)| psi);
            let got =
                orbital.nodaps_true_node_channel(fast_tt, request_ut1, constants::FLG_NONUT, &dpsi);
            if window.is_some_and(|(_, hi)| hi < request_ut1) {
                assert!(
                    matches!(got, Err(Error::EphemerisRange { requested_jd, body: Some(constants::TRUE_NODE), .. })
                    if requested_jd.to_bits() == request_ut1.to_bits())
                );
            } else {
                // Pinned v 3.2.1 on the matched partial fixture: slow geometry
                // at global TT, rather than the ~4.8-second-shifted fast TT.
                let node = got.expect("admitted request reaches slow lunar geometry");
                assert!((node[0] - 15.614141424545002).abs() <= 1e-9);
                assert_eq!(node[1].to_bits(), 0.0_f64.to_bits());
                for rate in &node[3..] {
                    assert_eq!(rate.to_bits(), 0.0_f64.to_bits());
                }
            }
        }
    }

    #[test]
    fn prefer_mean_low_byte_contract() {
        // (method, expected): mean when the low byte is zero or
        // NODBIT_MEAN is set; bit 8 never influences the choice.
        for (method, expected) in [
            (0, true),
            (1, true),
            (2, false),
            (3, true),
            (4, false),
            (6, false),
            (256, true),
            (257, true),
            (258, false),
            (1024, true),
            (-1, true),
            (-2, false),
        ] {
            assert_eq!(prefer_mean(method), expected, "method {method}");
        }
    }

    #[test]
    fn validate_orbital_object_branches() {
        // Valid planets pass.
        for body in [constants::MERCURY, constants::MARS, constants::PLUTO] {
            validate_orbital_object("get_orbital_elements", body, false).unwrap();
        }
        // Negative ids are rejected with the object-not-valid message.
        let err = validate_orbital_object("get_orbital_elements", -1, false).unwrap_err();
        assert_eq!(
            format!("{err}"),
            "get_orbital_elements: error in get_orbital_elements(): object -1 not valid"
        );
        // Lunar nodes/apogees and interpolated apsides are rejected.
        for body in [
            constants::MEAN_NODE,
            constants::TRUE_NODE,
            constants::MEAN_APOG,
            constants::OSCU_APOG,
            constants::INTP_APOG,
            constants::INTP_PERG,
        ] {
            validate_orbital_object("get_orbital_elements", body, false).unwrap_err();
        }
        // The Sun needs the sun_valid flag.
        validate_orbital_object("get_orbital_elements", constants::SUN, false).unwrap_err();
        validate_orbital_object("get_orbital_elements", constants::SUN, true).unwrap();
        // The 23-39 block carries its own message.
        let err = validate_orbital_object("get_orbital_elements", 25, false).unwrap_err();
        assert_eq!(
            format!("{err}"),
            "get_orbital_elements: illegal planet number 25."
        );
        // The sourceless block past the last fictitious body is rejected.
        validate_orbital_object("get_orbital_elements", constants::WALDEMATH + 1, false)
            .unwrap_err();
        // Harness category is input validation.
        assert_eq!(
            err.category(),
            crate::errors::ErrorCategory::InputValidation
        );
    }

    #[test]
    fn mean_elements_j2000_pins() {
        // Mercury semi-major axis is a stored constant at J2000.
        let mercury = mean_orbital_elements(constants::MERCURY, J2000_JD).unwrap();
        assert_eq!(mercury[0], 0.387098310);
        assert!((mercury[1] - 0.20563175).abs() < 1e-12);
        // Earth has no node: Omega and inclination are exactly zero.
        let earth = mean_orbital_elements(constants::EARTH, J2000_JD).unwrap();
        assert_eq!(earth[2], 0.0);
        assert_eq!(earth[3], 0.0);
        assert!((earth[0] - 1.000001018).abs() < 1e-12);
        // Pluto has no mean elements: callers fall back to osculating fits.
        assert!(mean_orbital_elements(constants::PLUTO, J2000_JD).is_none());
        // Non-finite epochs fail closed instead of propagating NaN.
        assert!(mean_orbital_elements(constants::MARS, f64::NAN).is_none());
    }

    #[test]
    fn semi_major_axis_circular_and_degenerate() {
        // Circular 1 AU orbit: v^2 = gm gives a = 1 exactly.
        let v = GM_SUN.sqrt();
        let a = osculating_semi_major_axis(1.0, 0.0, 0.0, 0.0, v, 0.0, GM_SUN);
        assert!((a - 1.0).abs() < 1e-15);
        // Degenerate zero state returns 0, never NaN.
        assert_eq!(
            osculating_semi_major_axis(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, GM_SUN),
            0.0
        );
    }

    #[test]
    fn conic_nodes_are_exact_reference_plane_intersections() {
        for inclination in [0.1, 0.8, 1.7, 2.8] {
            for node in [0.0, 0.7, 2.2, 5.9] {
                for argument in [0.0, 0.3, 1.2, 3.6, 6.0] {
                    let shape = OrbitShape {
                        a: 2.3,
                        e: 0.37,
                        incl: inclination,
                        node,
                        arg_peri: argument,
                        p: 2.3 * (1.0 - 0.37_f64.powi(2)),
                    };
                    let second_focus = focal_point_3d(shape).0;
                    for ascending in [true, false] {
                        let (point, radius) = orbit_node_3d(shape, ascending);
                        assert_eq!(point[2], 0.0);
                        let other_distance = point
                            .into_iter()
                            .zip(second_focus)
                            .map(|(coordinate, focus)| (coordinate - focus).powi(2))
                            .sum::<f64>()
                            .sqrt();
                        // An ellipse is the locus with constant sum of its
                        // distances from the two foci, independent of anomaly.
                        assert!(
                            (radius + other_distance - 2.0 * shape.a).abs()
                                <= 8.0 * f64::EPSILON * shape.a
                        );
                        let normal = [
                            node.sin() * inclination.sin(),
                            -node.cos() * inclination.sin(),
                            inclination.cos(),
                        ];
                        let plane_residual = point
                            .into_iter()
                            .zip(normal)
                            .map(|(coordinate, axis)| coordinate * axis)
                            .sum::<f64>();
                        assert!(plane_residual.abs() <= 4.0 * f64::EPSILON * radius);
                        let crossing = normal[0] * point[1] - normal[1] * point[0];
                        assert_eq!(crossing > 0.0, ascending);
                    }
                }
            }
        }
    }

    #[test]
    fn node_derivative_stays_on_owned_side_of_a_polynomial_join() {
        let epoch = J2000_JD;
        let (forward, h) = nodaps_derivative_stencil(epoch, epoch, epoch + 4.0).unwrap();
        assert_eq!(forward, [1.0, 2.0, 3.0, 4.0]);
        let (backward, _) = nodaps_derivative_stencil(epoch, epoch - 4.0, epoch).unwrap();
        assert_eq!(backward, [-1.0, -2.0, -3.0, -4.0]);
        for offsets in [forward, backward, [-2.0, -1.0, 1.0, 2.0]] {
            // The two adjacent records can have different constant/linear
            // terms: only this central record's derivative is meaningful.
            let trajectory = |t: f64| 2.0 + t / 8.0 + t * t + t.powi(4);
            let displacement = offsets.map(|n| trajectory(n * h) - trajectory(0.0));
            assert!((nodaps_stencil_derivative(displacement, offsets, h) - 0.125).abs() < 1e-13);
            assert_eq!(nodaps_stencil_derivative([0.0; 4], offsets, h), 0.0);
        }
    }

    #[test]
    fn node_derivative_respects_a_short_common_lunar_clock_domain() {
        let epoch = J2000_JD;
        let end = epoch + 5.0 / 86400.0;
        let (offsets, step) = nodaps_derivative_stencil(epoch, epoch, end).unwrap();
        assert!(step < NODAPS_SPEED_STEP_DAYS);
        assert!(offsets.into_iter().all(|offset| {
            let sample = epoch + offset * step;
            sample >= epoch && sample <= end && sample != epoch
        }));
    }

    #[test]
    fn spherical_to_cartesian_identity_and_radial_rate() {
        let cart = ecliptic_spherical_to_cartesian_state([0.0, 0.0, 1.0, 0.0, 0.0, 0.0]);
        assert_eq!(cart, [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        // Pure radial recession maps onto +x velocity alone.
        let cart = ecliptic_spherical_to_cartesian_state([0.0, 0.0, 2.0, 0.0, 0.0, 0.01]);
        assert!((cart[0] - 2.0).abs() < 1e-15);
        assert!((cart[3] - 0.01).abs() < 1e-15);
        assert_eq!(cart[4], 0.0);
        assert_eq!(cart[5], 0.0);
    }

    #[test]
    fn elements_of_circular_equatorial_orbit() {
        // Circular equatorial state at 1 AU: e = 0, i = 0, a = 1,
        // reserved tail slots stay zero.
        let v = GM_SUN.sqrt();
        let out = orbital_elements_from_ecliptic_state(
            [1.0, 0.0, 0.0, 0.0, v, 0.0],
            GM_SUN,
            constants::MARS,
            J2000_JD,
            None,
        );
        assert!((out[0] - 1.0).abs() < 1e-12);
        assert_eq!(out[1], 0.0);
        assert_eq!(out[2], 0.0);
        assert!((out[15] - 1.0).abs() < 1e-12);
        assert!((out[16] - 1.0).abs() < 1e-12);
        assert!(out[17..].iter().all(|s| *s == 0.0));
        // Mean motion matches the Gaussian constant in deg/day.
        assert!((out[11] - GAUSS_K.to_degrees()).abs() < 1e-12);
    }

    #[test]
    fn mean_shape_perihelion_aphelion_distances() {
        // Mercury J2000 mean shape: conic radius at nu=0/pi pins
        // a(1-e) and a(1+e).
        let el = mean_orbital_elements(constants::MERCURY, J2000_JD).unwrap();
        let shape = mean_shape_from_elements(el[0], el[1], el[2], el[3], el[4]);
        let (_, r_peri) = orbit_point_3d(shape, 0.0);
        let (_, r_aphe) = orbit_point_3d(shape, std::f64::consts::PI);
        assert!((r_peri - el[0] * (1.0 - el[1])).abs() < 1e-12);
        assert!((r_aphe - el[0] * (1.0 + el[1])).abs() < 1e-12);
        // Circular shape: node/apse longitudes coincide on the circle.
        let circ = mean_shape_from_elements(1.0, 0.0, 0.0, 0.0, 0.0);
        let (asc, r_asc) = orbit_point_3d(circ, -circ.arg_peri());
        assert!((r_asc - 1.0).abs() < 1e-15);
        assert!((asc[0] - 1.0).abs() + asc[1].abs() + asc[2].abs() < 1e-15);
    }

    #[test]
    fn osculating_shape_of_equatorial_circular_state() {
        let v = GM_SUN.sqrt();
        let shape = osculating_shape_from_ecliptic([1.0, 0.0, 0.0], [0.0, v, 0.0], GM_SUN);
        assert!((shape.semi_major_axis() - 1.0).abs() < 1e-12);
        assert!(shape.eccentricity() < 1e-10);
        assert_eq!(shape.inclination(), 0.0);
        assert_eq!(shape.node(), 0.0);
        assert_eq!(shape.arg_peri(), 0.0);
    }

    #[test]
    fn focal_point_second_focus_property() {
        // |focal| = 2 ae, carried in the aphelion slot under FOPOINT.
        let shape = mean_shape_from_elements(1.523679342, 0.09340062, 0.0, 0.0, 0.0);
        let (f, f_dist) = focal_point_3d(shape);
        assert!((f_dist - 2.0 * 1.523679342 * 0.09340062).abs() < 1e-12);
        let norm = (f[0] * f[0] + f[1] * f[1] + f[2] * f[2]).sqrt();
        assert!((norm - f_dist).abs() < 1e-15);
    }

    #[test]
    fn nodaps_ecliptic_height_is_independent_of_nutation_angles() {
        for epoch in [J2000_JD, 2460310.5, 2500000.25] {
            let vector = [0.3, -0.7, 0.31];
            let mean = NodapsFrame::for_date(epoch, None).unwrap();
            let height = mean.icrs_to_ecliptic(vector)[2];
            for nutation in [(0.0, 0.0), (0.01, -0.02), (-0.02, 0.03)] {
                let frame = NodapsFrame::for_date(epoch, Some(nutation)).unwrap();
                assert_eq!(frame.icrs_to_ecliptic(vector)[2], height);
                let round_trip = frame.ecliptic_to_icrs(frame.icrs_to_ecliptic(vector));
                assert!(
                    round_trip
                        .into_iter()
                        .zip(vector)
                        .all(|(a, b)| (a - b).abs() < 2e-14)
                );
            }
        }
    }

    #[test]
    fn frame_round_trip_and_j2000_x_axis() {
        let frame = NodapsFrame::for_date(2460311.0, None).unwrap();
        let v = [0.3, -0.7, 1.1];
        let back = frame.ecliptic_to_icrs(frame.icrs_to_ecliptic(v));
        // Two chained 3 x 3 rotations plus the obliquity pair accumulate
        // ordinary rounding (~20 flops at unit scale); measured 1.4 e-14.
        for (a, b) in back.into_iter().zip(v) {
            assert!((a - b).abs() < 2e-14);
        }
        // At J2000 the biased precession carries the ICRS/J2000
        // frame bias (~27 mas ≈ 1.3 e-7 rad observed), so the shared
        // x-axis is fixed only up to that scale — not identical.
        let j2k = NodapsFrame::for_date(J2000_JD, None).unwrap();
        let x = j2k.icrs_to_ecliptic([1.0, 0.0, 0.0]);
        assert!((x[0] - 1.0).abs() + x[1].abs() + x[2].abs() < 1e-6);
        // Non-finite epochs are refused.
        NodapsFrame::for_date(f64::NAN, None).unwrap_err();
    }

    #[test]
    fn j2000_vector_chain_is_identity_at_j2000() {
        let v = [0.2, 0.5, -0.8];
        let back = nodaps_vector_to_j2000(v, J2000_JD).unwrap();
        for (a, b) in back.into_iter().zip(v) {
            assert!((a - b).abs() < 1e-12);
        }
        let (lon, lat) = nodaps_dir_to_j2000(123.0, 4.5, J2000_JD).unwrap();
        assert!((lon - 123.0).abs() < 1e-9);
        assert!((lat - 4.5).abs() < 1e-9);
    }

    #[test]
    fn point_speeds_are_exact_for_a_cubic_trajectory_across_longitude_zero() {
        let sample = |t: f64| {
            (
                (0.005 + 0.5 * t + 0.1 * t.powi(2) + 0.01 * t.powi(3)).rem_euclid(360.0),
                0.0002 - 4e-6 * t + 0.1 * t.powi(2) - 0.2 * t.powi(3),
                1.0 + 0.01 * t + 0.02 * t.powi(2) + 0.03 * t.powi(3),
            )
        };
        let h = NODAPS_SPEED_STEP_DAYS;
        let result = nodaps_point_speeds(
            sample(0.0),
            sample(-h),
            sample(h),
            sample(-2.0 * h),
            sample(2.0 * h),
        );
        assert_eq!(&result[..3], &[0.005, 0.0002, 1.0]);
        assert!((result[3] - 0.5).abs() < 1e-11);
        assert!((result[4] - -4e-6).abs() < 1e-17);
        // Each cubic sample has at most eight rounded operations; the sum of
        // absolute derivative weights is 18/(12 h). Unlike the small latitude,
        // the distance has an O(1) offset whose sample rounding is amplified.
        let unit = f64::EPSILON / 2.0;
        let gamma = |n: f64| n * unit / (1.0 - n * unit);
        let extent = 2.0 * h;
        let distance_scale = 1.0 + 0.01 * extent + 0.02 * extent.powi(2) + 0.03 * extent.powi(3);
        let sample_bound = gamma(8.0) * distance_scale * 18.0 / (12.0 * h);
        let near = sample(h).2 - sample(-h).2;
        let far = sample(2.0 * h).2 - sample(-2.0 * h).2;
        let derivative_bound = gamma(3.0) * (8.0 * near.abs() + far.abs()) / (12.0 * h);
        assert!((result[5] - 0.01).abs() <= sample_bound + derivative_bound);
    }

    #[test]
    fn derivative_stencil_reduces_roundoff_in_a_nearly_stationary_latitude() {
        let h = NODAPS_SPEED_STEP_DAYS;
        let rate = -4e-6;
        let sample = |t: f64, error: f64| (10.0, 0.0002 + rate * t + error, 1.0);
        // A one-step position error must remain a small velocity error even
        // for rates much smaller than the coordinate's numerical resolution.
        let result = nodaps_point_speeds(
            sample(0.0, 0.0),
            sample(-h, -5e-15),
            sample(h, 5e-15),
            sample(-2.0 * h, -5e-15),
            sample(2.0 * h, 5e-15),
        );
        assert!((result[4] - rate).abs() < 4e-13);
        for jd in [-1_000_000.5_f64, 2451545.0, 8_000_000.5] {
            assert_eq!((jd + h) - (jd - h), 2.0 * h);
        }
    }

    #[test]
    fn leb_smooth_intervals_follow_actual_record_ownership() {
        let file = crate::leb::LebFile::open(&crate::test_data::core_leb_path()).unwrap();
        let index = file.index().unwrap();
        for body in [
            constants::MOON,
            constants::MERCURY,
            constants::EARTH,
            constants::TRUE_NODE,
        ] {
            let entry = index.find_body(body).unwrap();
            let join = entry.jd_start() + entry.interval_days();
            let left = entry.smooth_interval(join.next_down()).unwrap().unwrap();
            let right = entry.smooth_interval(join).unwrap().unwrap();
            assert!(left.1 < join, "body {body}: {left:?}");
            assert_eq!(right.0, join, "body {body}: {right:?}");
            assert!(right.1 > join);
        }
    }

    #[test]
    fn exact_nutation_rate_uses_the_central_record_at_a_join() {
        let file = crate::leb::LebFile::open(&crate::test_data::core_leb_path()).unwrap();
        let index = file.index().unwrap();
        let bytes = file.bytes().unwrap();
        let section = index.section(crate::leb::SECTION_NUTATION).unwrap();
        let header = crate::leb::NutationHeader::decode(bytes, section.offset() as usize).unwrap();
        let join = header.jd_start() + header.interval_days();
        let exact = index.nutation_rate(bytes, join).unwrap().unwrap();
        let h = 1.0 / 256.0;
        let current = index.eval_nutation(bytes, join).unwrap();
        let samples: Vec<_> = (1..=6)
            .map(|n| index.eval_nutation(bytes, join + n as f64 * h).unwrap())
            .collect();
        // Independent sixth-order endpoint interpolation tests a derivative
        // of the actual polynomial without stepping into its predecessor.
        for (component, rate) in [exact.0, exact.1].into_iter().enumerate() {
            let value = |pair: (f64, f64)| if component == 0 { pair.0 } else { pair.1 };
            let weights = [360.0, -450.0, 400.0, -225.0, 72.0, -10.0];
            let numerical = samples
                .iter()
                .zip(weights)
                .map(|(sample, weight)| (value(*sample) - value(current)) * weight)
                .sum::<f64>()
                / (60.0 * h);
            assert!(
                (rate - numerical).abs() < 2e-16,
                "component {component}: {rate} vs {numerical}"
            );
        }
    }

    #[test]
    fn mean_lunar_apsides_keep_the_model_radius_and_zero_radial_derivative() {
        let runtime = crate::config::Runtime::from_context(
            crate::context::EphemerisContext::new()
                .with_mode(crate::context::ProviderMode::SealedLeb),
        );
        runtime.set_leb_file(Some(crate::test_data::core_leb_path()));
        runtime.set_iers_delta_t_enabled(Some(false));
        let session = runtime.session().unwrap();
        for epoch in [2451545.0, 2460310.5, 2460310.75] {
            let points = session
                .nod_aps(
                    epoch,
                    constants::MOON,
                    constants::NODBIT_MEAN,
                    constants::FLG_SWIEPH | constants::FLG_SPEED,
                )
                .unwrap();
            let points = [
                points.ascending,
                points.descending,
                points.perihelion,
                points.aphelion,
            ];
            assert_eq!(
                points[2][2],
                2.0 * MEAN_LUNAR_DISTANCE_AU - crate::lunar::MEAN_APOGEE_DISTANCE_AU
            );
            assert_eq!(points[3][2], crate::lunar::MEAN_APOGEE_DISTANCE_AU);
            assert_eq!(points[2][5], 0.0);
            assert_eq!(points[3][5], 0.0);
            // Both conic nodes obey 1/r_asc + 1/r_desc = 2/p, so its
            // differentiated identity is zero at every phase, including
            // a discontinuity in the independently stored nutation curve.
            let reciprocal_rate =
                points[0][5] / points[0][2].powi(2) + points[1][5] / points[1][2].powi(2);
            assert!(reciprocal_rate.abs() < 2e-14);
            assert_eq!(points[2][3], points[3][3]);
            assert_eq!(points[2][4], -points[3][4]);
        }
    }

    #[test]
    fn node_latitude_rate_matches_a_higher_order_position_derivative() {
        let runtime = crate::config::Runtime::from_context(
            crate::context::EphemerisContext::new()
                .with_mode(crate::context::ProviderMode::SealedLeb),
        );
        runtime.set_leb_file(Some(crate::test_data::core_leb_path()));
        runtime.set_iers_delta_t_enabled(Some(false));
        runtime.set_delta_t_userdef(Some(64.0 / 86400.0));
        let session = runtime.session().unwrap();
        let ut_tt = crate::time::default_timescale_tt_from_ut1(2451545.0).unwrap();
        for (epoch, body, method) in [
            (2451545.0, constants::MERCURY, constants::NODBIT_MEAN),
            (2451545.0, constants::MERCURY, constants::NODBIT_OSCU),
            (ut_tt, constants::MARS, constants::NODBIT_OSCU),
        ] {
            let rate = session
                .nod_aps(
                    epoch,
                    body,
                    method,
                    constants::FLG_SWIEPH | constants::FLG_SPEED,
                )
                .unwrap()
                .descending[4];
            // Independent sixth-order differentiation of the public positions
            // at a different spacing tests velocity = d(position)/d(TT).
            // No previously stored output or ephemeris-derived rate is used.
            let h = 1.0 / 32.0;
            let latitude = |offset: f64| {
                session
                    .nod_aps(epoch + offset, body, method, constants::FLG_SWIEPH)
                    .unwrap()
                    .descending[1]
            };
            let difference = |n: f64| latitude(n * h) - latitude(-n * h);
            let derivative =
                (45.0 * difference(1.0) - 9.0 * difference(2.0) + difference(3.0)) / (60.0 * h);
            assert!(
                (rate - derivative).abs() < 2e-12,
                "body {body}, method {method}: rate {rate}, position derivative {derivative}"
            );
        }
    }

    #[test]
    fn intp_sentinel_shape() {
        let plain = intp_undefined_points(false);
        for pt in plain {
            assert!(pt[0..3].iter().all(|s| s.is_nan()));
            assert_eq!(&pt[3..], &[0.0, 0.0, 0.0]);
        }
        let speed = intp_undefined_points(true);
        for pt in speed {
            assert!(pt.iter().all(|s| s.is_nan()));
        }
    }

    #[test]
    fn moon_composition_geocentric_mean_core() {
        // Geocentric mean composition: node longitude preserved,
        // perigee opposite the apogee, apogee distance preserved.
        let pts = moon_nodaps_points(
            100.0, 5.0, 200.0, -3.0, 0.0027, 0.0024, false, false, false, [0.0; 3], [0.0; 3],
            false, 2460311.0, None,
        )
        .unwrap();
        assert!((pts[0][0] - 100.0).abs() < 1e-9);
        assert!((pts[0][1] - 5.0).abs() < 1e-9);
        assert!((pts[2][0] - 20.0).abs() < 1e-9);
        assert!((pts[2][1] - 3.0).abs() < 1e-9);
        assert!((pts[2][2] - 0.0024).abs() < 1e-15);
        assert!((pts[3][2] - 0.0027).abs() < 1e-15);
        // Focal substitution carries apo - peri in the aphelion slot.
        let focal = moon_nodaps_points(
            100.0, 5.0, 200.0, -3.0, 0.0027, 0.0024, true, false, false, [0.0; 3], [0.0; 3], false,
            2460311.0, None,
        )
        .unwrap();
        assert!((focal[3][2] - 0.0003).abs() < 1e-15);
        // Non-finite and non-positive inputs are refused.
        moon_nodaps_points(
            f64::NAN,
            5.0,
            200.0,
            -3.0,
            0.0027,
            0.0024,
            false,
            false,
            false,
            [0.0; 3],
            [0.0; 3],
            false,
            2460311.0,
            None,
        )
        .unwrap_err();
        moon_nodaps_points(
            100.0, 5.0, 200.0, -3.0, 0.0, 0.0024, false, false, false, [0.0; 3], [0.0; 3], false,
            2460311.0, None,
        )
        .unwrap_err();
    }

    /// Serving-level Sun rejection: both public element
    /// entries reject the Sun with the model `sun_valid = false`
    /// message, never the zero tuple. Validation runs before any
    /// provider read, so an unreachable stub provider suffices.
    #[test]
    fn sun_elements_refused_like_the_release() {
        struct NoProvider;
        impl crate::positions::StateProvider for NoProvider {
            fn provider_name(&self) -> &'static str {
                "test/stub"
            }
            fn has_body(&self, _body_id: i32) -> bool {
                unreachable!("validation precedes the provider read")
            }
            fn body_coverage(&self, _body_id: i32) -> Option<(f64, f64)> {
                unreachable!("validation precedes the provider read")
            }
            fn serving_coord_type(&self, _body_id: i32) -> Option<u32> {
                unreachable!("validation precedes the provider read")
            }
            fn eval_raw(&self, _body_id: i32, _jd: f64) -> Result<(Vec<f64>, Vec<f64>)> {
                unreachable!("validation precedes the provider read")
            }
            fn eval_raw_two_part(
                &self,
                _body_id: i32,
                _epoch: crate::types::TwoPartEpoch,
            ) -> Result<(Vec<f64>, Vec<f64>)> {
                unreachable!("validation precedes the provider read")
            }
            fn serving_entry_snapshot(
                &self,
                _body_id: i32,
            ) -> Result<crate::positions::EntrySnapshot> {
                unreachable!("validation precedes the provider read")
            }
        }
        let provider = NoProvider;
        let nutation_at = |_: f64| -> Result<(f64, f64)> {
            unreachable!("validation precedes the nutation read")
        };
        let serving = crate::positions::SealedServing::new(&provider, &nutation_at);
        let deltat_at =
            |_: f64| -> Result<f64> { unreachable!("validation precedes the Delta-T read") };
        let bound = SealedOrbitalServing::new(&serving, &deltat_at);
        let err = bound
            .orbital_elements_ut(constants::SUN, 2460311.0, 0)
            .unwrap_err();
        assert_eq!(
            format!("{err}"),
            "get_orbital_elements_ut: error in get_orbital_elements(): object 0 not valid"
        );
        let err = bound
            .orbital_elements_tt(constants::SUN, 2460311.0, 0)
            .unwrap_err();
        assert_eq!(
            format!("{err}"),
            "get_orbital_elements: error in get_orbital_elements(): object 0 not valid"
        );
    }

    #[test]
    fn invalid_orbital_object_error_contract() {
        let err = Error::InvalidOrbitalObject {
            entry: "nod_aps",
            detail: "nodes/apsides for body id 59 are not implemented.".to_string(),
        };
        assert_eq!(
            format!("{err}"),
            "nod_aps: nodes/apsides for body id 59 are not implemented."
        );
    }
}

#[cfg(test)]
mod body02_tests {
    use super::*;

    fn circle(a_au: f64) -> OrbitEllipse {
        OrbitEllipse::new(a_au, 0.0, 0.0, 0.0, 0.0).unwrap()
    }

    fn elements_with(a: f64, e: f64) -> [f64; 50] {
        let mut el = [0.0f64; 50];
        el[0] = a;
        el[1] = e;
        el[15] = a * (1.0 - e);
        el[16] = a * (1.0 + e);
        el
    }

    #[test]
    fn ellipse_position_circular_orbit() {
        // Circular equatorial unit orbit: cardinal anomalies land on the
        // axes and the radius is identically 1.
        let el = circle(1.0);
        for (anom, want) in [
            (0.0, [1.0, 0.0, 0.0]),
            (std::f64::consts::FRAC_PI_2, [0.0, 1.0, 0.0]),
            (std::f64::consts::PI, [-1.0, 0.0, 0.0]),
        ] {
            let p = ellipse_position_au(el, anom);
            for k in 0..3 {
                assert!((p[k] - want[k]).abs() < 1e-15, "k={k} p={p:?}");
            }
        }
        // Eccentric orbit: perihelion/aphelion distances on the apsidal line.
        let el = OrbitEllipse::new(1.5, 0.2, 0.0, 0.0, 0.0).unwrap();
        let peri = ellipse_position_au(el, 0.0);
        let apo = ellipse_position_au(el, std::f64::consts::PI);
        assert!((peri[0] - 1.2).abs() < 1e-15);
        assert!((apo[0] + 1.8).abs() < 1e-15);
    }

    #[test]
    fn helio_extremes_read_peri_apo_slots() {
        let (max, min) = helio_distance_extremes(&elements_with(1.5, 0.2)).unwrap();
        assert!((max - 1.8).abs() < 1e-15, "max={max}");
        assert!((min - 1.2).abs() < 1e-15, "min={min}");
        // Degenerate semi-major axis: the model zero triple.
        assert_eq!(helio_distance_extremes(&elements_with(0.0, 0.0)), None);
        assert_eq!(helio_distance_extremes(&elements_with(-1.0, 0.5)), None);
    }

    #[test]
    fn coaxial_circles_min_max() {
        // Two circular equatorial orbits: min at conjunction (0.5 AU),
        // max at opposition (2.5 AU), refined well under 1 e-8 AU.
        let (max, min) = two_ellipse_min_max(circle(1.0), circle(1.5));
        assert!((min - 0.5).abs() < 1e-9, "min={min}");
        assert!((max - 2.5).abs() < 1e-9, "max={max}");
    }

    #[test]
    fn geocentric_extremes_degenerate_contract() {
        // Invalid or non-elliptical osculating orbit: the model zero triple.
        assert_eq!(
            geocentric_two_ellipse_extremes(&elements_with(0.0, 0.0), &elements_with(1.0, 0.0))
                .unwrap(),
            None
        );
        assert_eq!(
            geocentric_two_ellipse_extremes(&elements_with(1.5, 1.0), &elements_with(1.0, 0.0))
                .unwrap(),
            None
        );
        assert_eq!(
            geocentric_two_ellipse_extremes(&elements_with(1.5, 1.2), &elements_with(1.0, 0.0))
                .unwrap(),
            None
        );
        // Valid pair agrees with the direct solver.
        let direct = two_ellipse_min_max(circle(1.0), circle(1.5));
        assert_eq!(
            geocentric_two_ellipse_extremes(&elements_with(1.0, 0.0), &elements_with(1.5, 0.0))
                .unwrap(),
            Some(direct)
        );
    }

    #[test]
    fn sun_earth_emb_composition() {
        let emb = elements_with(1.0, 0.0167);
        assert_eq!(sun_extremes_from_emb_elements(&emb), (1.0167, 0.9833));
        assert_eq!(earth_extremes_from_emb_elements(&emb), (2.0, 0.0));
    }

    #[test]
    fn emb_channel_composition_matches_release_order() {
        // : earth + (moon - earth) * w
        // per component with w = 1 / (81.3005691 + 1).
        let w = emb_composition_weight();
        assert!((w - 1.0 / 82.3005691).abs() == 0.0);
        let earth = ([1.0, 0.0, 0.0], [0.0, 0.01, 0.0]);
        let moon = ([1.00257, 0.0, 0.0], [0.0, 0.0105, 0.0]);
        let (p, v) = emb_icrs_state_from_channels(earth, moon);
        assert!((p[0] - (1.0 + 0.00257 * w)).abs() < 1e-16);
        assert_eq!(p[1], 0.0);
        assert!((v[1] - (0.01 + 0.0005 * w)).abs() < 1e-18);
        // Coincident channels stay fixed (the Earth arm minimum).
        let (p0, _) = emb_icrs_state_from_channels(earth, earth);
        assert_eq!(p0, earth.0);
    }

    #[test]
    fn orbit_ellipse_from_elements_degrees_to_radians() {
        let mut el = [0.0f64; 50];
        el[0] = 2.0;
        el[1] = 0.1;
        el[2] = 90.0;
        el[3] = 180.0;
        el[4] = 45.0;
        let e = orbit_ellipse_from_elements(&el).unwrap();
        assert!((e.a_au() - 2.0).abs() == 0.0);
        assert!((e.inc_rad() - std::f64::consts::FRAC_PI_2).abs() < 1e-16);
        assert!((e.node_rad() - std::f64::consts::PI).abs() < 1e-16);
        assert!((e.arg_peri_rad() - std::f64::consts::FRAC_PI_4).abs() < 1e-16);
        // Non-finite values are refused, never searched.
        el[0] = f64::NAN;
        orbit_ellipse_from_elements(&el).unwrap_err();
    }

    #[test]
    fn flag_helpers_match_release_heads() {
        // FLG_HELCTR is a no-op for the Sun and the Moon only.
        assert_eq!(
            strip_helio_for_sun_moon(constants::SUN, constants::FLG_HELCTR),
            0
        );
        assert_eq!(
            strip_helio_for_sun_moon(constants::MOON, constants::FLG_HELCTR),
            0
        );
        assert_eq!(
            strip_helio_for_sun_moon(constants::MARS, constants::FLG_HELCTR),
            constants::FLG_HELCTR
        );
        // Interior BARYCTR degrades to HELCTR; trans-jovian survives.
        let bary = constants::FLG_BARYCTR;
        let helio = constants::FLG_HELCTR;
        assert_eq!(
            normalize_barycentric_flag(constants::MARS, bary, 1.52),
            helio
        );
        assert_eq!(
            normalize_barycentric_flag(constants::JUPITER, bary, 5.20),
            helio
        );
        assert_eq!(
            normalize_barycentric_flag(constants::SATURN, bary, 9.58),
            bary
        );
        assert_eq!(normalize_barycentric_flag(constants::MARS, 0, 1.52), 0);
        // Geocentric fictitious constructions.
        assert!(is_fict_geocentric_body(constants::WHITE_MOON));
        assert!(is_fict_geocentric_body(constants::WALDEMATH));
        assert!(!is_fict_geocentric_body(constants::MARS));
        // Validation: the Sun is valid here, 23-39 and negatives are not.
        validate_orbital_object("orbit_max_min_true_distance", constants::SUN, true).unwrap();
        validate_orbital_object("orbit_max_min_true_distance", 25, true).unwrap_err();
        validate_orbital_object("orbit_max_min_true_distance", -2, true).unwrap_err();
        validate_orbital_object(
            "orbit_max_min_true_distance",
            constants::WALDEMATH + 1,
            true,
        )
        .unwrap_err();
    }
}

#[cfg(test)]
mod native_nodaps_review_tests {
    use super::*;
    use std::cell::RefCell;

    struct NoGeometry;
    impl StateProvider for NoGeometry {
        fn serving_entry_snapshot(&self, _: i32) -> Result<crate::positions::EntrySnapshot> {
            panic!("no stored entry")
        }
        fn provider_name(&self) -> &'static str {
            "test/native-nodaps-site"
        }
        fn eval_raw_two_part(
            &self,
            _: i32,
            _: crate::types::TwoPartEpoch,
        ) -> Result<(Vec<f64>, Vec<f64>)> {
            panic!("no stored split read")
        }
        fn has_body(&self, _: i32) -> bool {
            panic!("no stored support")
        }
        fn body_coverage(&self, _: i32) -> Option<(f64, f64)> {
            panic!("no stored domain")
        }
        fn serving_coord_type(&self, _: i32) -> Option<u32> {
            panic!("no stored frame")
        }
        fn eval_raw(&self, _: i32, _: f64) -> Result<(Vec<f64>, Vec<f64>)> {
            panic!("no stored read")
        }
    }

    #[test]
    fn native_site_undefined_points_acquire_conditional_states_per_point() {
        let provider = NoGeometry;
        let nut = |_: f64| -> Result<(f64, f64)> { panic!("NONUT") };
        let serving = SealedServing::new(&provider, &nut);
        let delta = |_: f64| -> Result<f64> { panic!("site observation has no model clock") };
        let tt = 2460400.25;
        for flags in [
            constants::FLG_NONUT,
            constants::FLG_NONUT | constants::FLG_SPEED,
            constants::FLG_NONUT | constants::FLG_TRUEPOS,
            constants::FLG_NONUT | constants::FLG_NOABERR,
        ] {
            let reads = RefCell::new(Vec::new());
            let site_at = |epoch: f64, _: (f64, f64, f64)| {
                reads.borrow_mut().push(epoch.to_bits());
                Ok(([0.00003, 0.00002, 0.00001], [0.0001, 0.0002, 0.0]))
            };
            let engine =
                SealedOrbitalServing::new(&serving, &delta).with_nodaps_site_state(&site_at);
            let points = engine
                .nod_aps_topo_tt(
                    constants::INTP_APOG,
                    tt,
                    tt,
                    1,
                    flags | constants::FLG_TOPOCTR,
                    (13.405, 52.52, 34.0),
                )
                .unwrap();
            let expected_point_reads =
                if flags & (constants::FLG_TRUEPOS | constants::FLG_NOABERR) == 0 {
                    vec![
                        tt.to_bits(),
                        tt.to_bits(),
                        (tt - 0.0001).to_bits(),
                        (tt + 0.0001).to_bits(),
                    ]
                } else {
                    vec![tt.to_bits()]
                };
            assert_eq!(*reads.borrow(), expected_point_reads.repeat(4));
            for point in points {
                assert!(point[0].is_nan() && point[2].is_nan());
                assert!(point[1].is_nan()); // Undefined geometry remains undefined.
                assert_eq!(point[3].to_bits(), 0.0_f64.to_bits());
                assert_eq!(point[4].to_bits(), 0.0_f64.to_bits());
                if flags & constants::FLG_SPEED != 0 {
                    assert!(point[5].is_nan());
                } else {
                    assert_eq!(point[5].to_bits(), 0.0_f64.to_bits());
                }
            }
        }
    }

    #[test]
    fn native_site_zero_and_centered_points_do_not_acquire_support() {
        let provider = NoGeometry;
        let nut = |_: f64| -> Result<(f64, f64)> { panic!("no frame needed") };
        let serving = SealedServing::new(&provider, &nut);
        let delta = |_: f64| -> Result<f64> { panic!("no model clock") };
        let site_at = |_: f64, _: (f64, f64, f64)| -> Result<([f64; 3], [f64; 3])> {
            panic!("zero/centered guard")
        };
        let engine = SealedOrbitalServing::new(&serving, &delta).with_nodaps_site_state(&site_at);
        let zero = [[0.0; 6]; 4];
        assert_eq!(
            engine
                .apply_nodaps_topocentric(
                    zero,
                    2451545.0,
                    0.0,
                    constants::FLG_TOPOCTR,
                    (0.0, 0.0, 0.0)
                )
                .unwrap(),
            zero
        );
        for center in [constants::FLG_HELCTR, constants::FLG_BARYCTR] {
            let points = [[20.0, 1.0, 1.0, 0.0, 0.0, 0.0]; 4];
            assert_eq!(
                engine
                    .apply_nodaps_topocentric(
                        points,
                        2451545.0,
                        0.0,
                        constants::FLG_TOPOCTR | center,
                        (0.0, 0.0, 0.0)
                    )
                    .unwrap(),
                points
            );
        }
    }

    #[test]
    fn apparent_site_side_domain_refuses_even_without_speed() {
        let provider = NoGeometry;
        let nut = |_: f64| -> Result<(f64, f64)> { panic!("NONUT") };
        let serving = SealedServing::new(&provider, &nut);
        let delta = |_: f64| -> Result<f64> { panic!("no model clock") };
        let tt = 2460400.25;
        let site_at = |epoch: f64, _: (f64, f64, f64)| {
            if epoch < tt {
                Err(Error::SpkOutOfRange {
                    target: 399,
                    et: -1.0,
                    start: 0.0,
                    end: 1.0,
                })
            } else {
                Ok(([0.00003, 0.00002, 0.00001], [0.0001, 0.0002, 0.0]))
            }
        };
        let engine = SealedOrbitalServing::new(&serving, &delta).with_nodaps_site_state(&site_at);
        let flags = constants::FLG_TOPOCTR | constants::FLG_NONUT;
        assert!(matches!(
            engine.nod_aps_topo_tt(constants::INTP_APOG, tt, tt, 1, flags, (0.0, 0.0, 0.0)),
            Err(Error::SpkOutOfRange {
                target: 399,
                et: -1.0,
                ..
            })
        ));
        for suppress in [constants::FLG_TRUEPOS, constants::FLG_NOABERR] {
            assert!(
                engine
                    .nod_aps_topo_tt(
                        constants::INTP_APOG,
                        tt,
                        tt,
                        1,
                        flags | suppress,
                        (0.0, 0.0, 0.0)
                    )
                    .is_ok()
            );
        }
    }

    struct RecordingAnchor(RefCell<Vec<(u64, bool, bool)>>);
    impl crate::sidereal::DynamicAnchorProvider for RecordingAnchor {
        fn anchor_longitude_deg(
            &self,
            _: crate::sidereal::DynamicAnchor,
            tt: f64,
            noaberr: bool,
            nogdefl: bool,
        ) -> Result<f64> {
            self.0.borrow_mut().push((tt.to_bits(), noaberr, nogdefl));
            Ok(180.0 + (tt - 2451545.0) * 0.001)
        }
    }

    #[test]
    fn native_body_and_site_ayanamsha_use_global_observation_tt_and_flag_toggles() {
        let provider = NoGeometry;
        let nut = |_: f64| -> Result<(f64, f64)> { panic!("no stored nutation") };
        let serving = SealedServing::new(&provider, &nut);
        let delta = |_: f64| Ok(0.25); // Distinguishes caller model from global Time.
        let original_ut = 2460400.25;
        let original_tt = crate::time::default_timescale_tt_from_ut1(original_ut).unwrap();
        let ut_at = |tt: f64| {
            if tt == original_tt {
                Ok(original_ut)
            } else {
                crate::time::default_timescale_ut1_from_tt(tt)
            }
        };
        let anchor = RecordingAnchor(RefCell::new(Vec::new()));
        let mut selection = crate::sidereal::SiderealSelection::unset();
        assert!(
            selection
                .set_sidereal_mode(constants::SIDM_TRUE_CITRA, 0.0, 0.0)
                .is_empty()
        );
        let aya = SealedAyanamshaServing::new(selection, &delta, Some(&anchor));
        let mut engine = SealedOrbitalServing::new(&serving, &delta).with_ayanamsha(aya);
        engine.nodaps_ut1_at = Some(&ut_at); // Exercise only the private clock tail.
        for flags in [
            0,
            constants::FLG_NOABERR,
            constants::FLG_NOGDEFL,
            constants::FLG_TRUEPOS,
            constants::FLG_HELCTR,
            constants::FLG_BARYCTR,
        ] {
            let stripped = (flags & !constants::FLG_J2000) | constants::FLG_NONUT;
            let (noaberr, nogdefl, _) = crate::sidereal::resolve_ayanamsha_flags(stripped);
            for tt in [original_tt, original_tt - 0.001, original_tt + 0.001] {
                let ut = ut_at(tt).unwrap();
                let expected_tt = crate::time::default_timescale_tt_from_ut1(ut).unwrap();
                engine.nodaps_mean_ayanamsha_deg(tt, flags).unwrap();
                assert_eq!(
                    anchor.0.borrow().last().copied(),
                    Some((expected_tt.to_bits(), noaberr, nogdefl))
                );
            }
            let site_ut = crate::time::default_timescale_ut1_from_tt(original_tt).unwrap();
            let expected_site_tt = crate::time::default_timescale_tt_from_ut1(site_ut).unwrap();
            engine
                .nodaps_true_ayanamsha_deg(original_tt, flags)
                .unwrap();
            assert_eq!(
                anchor.0.borrow().last().copied(),
                Some((expected_site_tt.to_bits(), noaberr, nogdefl))
            );
        }
    }

    #[test]
    fn native_ayanamsha_retains_configured_user_ut_and_valens_defining_epochs() {
        let provider = NoGeometry;
        let nut = |_: f64| -> Result<(f64, f64)> { panic!("no stored nutation") };
        let serving = SealedServing::new(&provider, &nut);
        let delta = |_: f64| Ok(0.25);
        let ut = 2460400.25;
        let global_tt = crate::time::default_timescale_tt_from_ut1(ut).unwrap();
        let clock = |_: f64| Ok(ut);
        let mut engine = SealedOrbitalServing::new(&serving, &delta);
        engine.nodaps_ut1_at = Some(&clock);
        for mode in [
            constants::SIDM_USER | constants::SIDBIT_USER_UT,
            constants::SIDM_VALENS_MOON,
        ] {
            let mut selection = crate::sidereal::SiderealSelection::unset();
            assert!(
                selection
                    .set_sidereal_mode(mode, 2451545.0, 23.0)
                    .is_empty()
            );
            let aya = SealedAyanamshaServing::new(selection, &delta, None);
            let actual = engine
                .nodaps_ayanamsha_at_ut1(&aya, ut, constants::FLG_NONUT)
                .unwrap();
            let expected = aya.mean_tt_with(global_tt, false, false).unwrap().0;
            assert_eq!(actual.to_bits(), expected.to_bits());
            let zero = |_: f64| Ok(0.0);
            let other = SealedAyanamshaServing::new(selection, &zero, None)
                .mean_tt_with(global_tt, false, false)
                .unwrap()
                .0;
            assert_ne!(actual.to_bits(), other.to_bits());
        }
    }
}
